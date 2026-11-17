package no.uio.imaging;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import ij.ImagePlus;
import ij.ImageStack;
import ij.process.FloatProcessor;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import org.junit.Test;

/**
 * Tests on synthetic frames: isolated Gaussian emitters (PSF sigma 130 nm) on a constant background with Poisson shot
 * noise and Gaussian read noise. Positions are compared with the ground truth.
 */
public class SMLM_LocalizationTest {

    static final double PSF_NM = 130.0;

    /** Ground-truth emitters on a jittered grid, in nm, for a field of fieldNm x fieldNm. */
    static double[][] emitters(double fieldNm, double spacingNm, long seed) {
        Random rnd = new Random(seed);
        List<double[]> e = new ArrayList<double[]>();
        for (double y = 1500; y < fieldNm - 1500; y += spacingNm)
            for (double x = 1500; x < fieldNm - 1500; x += spacingNm)
                e.add(new double[]{x + (rnd.nextDouble() - 0.5) * 400, y + (rnd.nextDouble() - 0.5) * 400});
        return e.toArray(new double[0][]);
    }

    /** One frame with the emitters (positions in nm, top-left corner = 0 nm) imaged at pixelNm. */
    static float[] frame(double[][] em, double pixelNm, int W, int H, double photons, double bg, double readSd,
                         double dxNm, double dyNm, Random rnd) {
        float[] img = new float[W * H];
        double s = PSF_NM / pixelNm;
        for (int y = 0; y < H; y++) {
            for (int x = 0; x < W; x++) {
                double v = bg;
                for (double[] e : em) {
                    double ex = (e[0] + dxNm) / pixelNm - 0.5, ey = (e[1] + dyNm) / pixelNm - 0.5;
                    if (Math.abs(ex - x) > 6 * s || Math.abs(ey - y) > 6 * s) continue;
                    // pixel-integrated Gaussian
                    v += photons * 0.25 * (erf((x + 0.5 - ex) / (Math.sqrt(2) * s)) - erf((x - 0.5 - ex) / (Math.sqrt(2) * s)))
                                        * (erf((y + 0.5 - ey) / (Math.sqrt(2) * s)) - erf((y - 0.5 - ey) / (Math.sqrt(2) * s)));
                }
                img[y * W + x] = (float) (poisson(v, rnd) + 100 + readSd * rnd.nextGaussian());
            }
        }
        return img;
    }

    static double poisson(double lam, Random rnd) {
        if (lam > 50) return Math.max(0, Math.round(lam + Math.sqrt(lam) * rnd.nextGaussian()));
        double L = Math.exp(-lam), p = 1; int k = 0;
        do { k++; p *= rnd.nextDouble(); } while (p > L);
        return k - 1;
    }

    static double erf(double x) {   // Abramowitz-Stegun 7.1.26, |error| < 1.5e-7
        double t = 1 / (1 + 0.3275911 * Math.abs(x));
        double y = 1 - (((((1.061405429 * t - 1.453152027) * t) + 1.421413741) * t - 0.284496736) * t + 0.254829592) * t * Math.exp(-x * x);
        return x >= 0 ? y : -y;
    }

    /** Matches localisations to the truth (nearest within 150 nm); returns {recall, rmse_nm, falsePositives}. */
    static double[] score(List<SMLM_Localization.Loc> locs, double[][] em, double pixelNm, double dxNm, double dyNm) {
        int found = 0; double se = 0; int fp = 0;
        boolean[] used = new boolean[em.length];
        for (SMLM_Localization.Loc l : locs) {
            double x = (l.x + 0.5) * pixelNm - dxNm, y = (l.y + 0.5) * pixelNm - dyNm;
            int best = -1; double bd = 150;
            for (int i = 0; i < em.length; i++) {
                double d = Math.hypot(em[i][0] - x, em[i][1] - y);
                if (d < bd) { bd = d; best = i; }
            }
            if (best < 0 || used[best]) { fp++; continue; }
            used[best] = true; found++; se += bd * bd;
        }
        return new double[]{found / (double) em.length, Math.sqrt(se / Math.max(1, found)), fp};
    }

    static SMLM_Localization.Params params(double pixelNm) {
        SMLM_Localization.Params p = new SMLM_Localization.Params();
        p.pixelNm = pixelNm;
        return p;
    }

    @Test
    public void pixelConversionMatchesReference() {
        SMLM_Localization.Params p = params(107);
        assertEquals(4, p.fitRadiusPx());        // 9 x 9 window at 107 nm/px
        assertEquals(3, p.minSeparationPx());    // 7 x 7 neighbourhood
        p = params(160);
        assertEquals(3, p.fitRadiusPx());        // 7 x 7 window at 160 nm/px
        assertEquals(2, p.minSeparationPx());    // 5 x 5 neighbourhood
        assertEquals(0.66875, p.dogSigma1Px(), 1e-9);
    }

    @Test
    public void gaussianFilterMatchesDirectConvolution() {
        int W = 23, H = 17; Random rnd = new Random(3);
        float[] img = new float[W * H];
        for (int i = 0; i < img.length; i++) img[i] = (float) rnd.nextGaussian();
        double s = 1.7; int rad = (int) (4 * s + 0.5);
        float[] f = SMLM_Localization.gaussianFilter(img, W, H, s);
        double norm = 0;
        for (int i = -rad; i <= rad; i++) norm += Math.exp(-0.5 * i * i / (s * s));
        int x = 0, y = 5;   // a border pixel exercises the mirror boundary
        double acc = 0;
        for (int j = -rad; j <= rad; j++) for (int i = -rad; i <= rad; i++)
            acc += Math.exp(-0.5 * (i * i + j * j) / (s * s)) / (norm * norm)
                 * img[SMLM_Localization.reflect(y + j, H) * W + SMLM_Localization.reflect(x + i, W)];
        assertEquals(acc, f[y * W + x], 1e-5);
    }

    @Test
    public void localisesSyntheticEmittersAt107And160nm() {
        double field = 16000;
        double[][] em = emitters(field, 1500, 11);
        for (double pix : new double[]{107, 160}) {
            int W = (int) Math.ceil(field / pix), H = W;
            Random rnd = new Random(5);
            List<SMLM_Localization.Loc> locs = new ArrayList<SMLM_Localization.Loc>();
            SMLM_Localization.Params p = params(pix);
            for (int f = 0; f < 5; f++)
                SMLM_Localization.localizeFrame(frame(em, pix, W, H, 1500, 20, 2, 0, 0, rnd), W, H, f, p, locs);
            List<SMLM_Localization.Loc> f0 = new ArrayList<SMLM_Localization.Loc>();
            for (SMLM_Localization.Loc l : locs) if (l.frame == 0) f0.add(l);
            double[] sc = score(f0, em, pix, 0, 0);
            double medSigma = 0; for (SMLM_Localization.Loc l : f0) medSigma += l.sigma * pix / f0.size();
            System.out.printf(java.util.Locale.ROOT, "pixel %.0f nm: %d emitters, recall %.3f, RMSE %.1f nm, FP %d, mean sigma %.1f nm%n",
                    pix, em.length, sc[0], sc[1], (int) sc[2], medSigma);
            assertTrue("recall " + sc[0], sc[0] > 0.97);
            assertTrue("rmse " + sc[1], sc[1] < 15);
            assertTrue("false positives " + sc[2], sc[2] <= 2);
            assertEquals(PSF_NM, medSigma, 12);
        }
    }

    @Test
    public void driftCorrectionRemovesImposedDriftWithCorrectSign() {
        double pix = 107, field = 12000;
        int W = (int) Math.ceil(field / pix), H = W, nF = 100;
        double[][] em = emitters(field, 700, 21);
        Random rnd = new Random(7);
        ImageStack st = new ImageStack(W, H);
        for (int f = 0; f < nF; f++) {
            // linear drift of +300 nm in x and -200 nm in y over the series; each frame shows a random 30 % of emitters
            List<double[]> on = new ArrayList<double[]>();
            for (double[] e : em) if (rnd.nextDouble() < 0.3) on.add(e);
            double t = f / (double) (nF - 1);
            st.addSlice(new FloatProcessor(W, H, frame(on.toArray(new double[0][]), pix, W, H, 1500, 20, 2, 300 * t, -200 * t, rnd), null));
        }
        ImagePlus imp = new ImagePlus("drift", st);
        SMLM_Localization.Params p = params(pix);
        p.driftCorrection = true; p.driftSegment = 10;
        List<SMLM_Localization.Loc> locs = SMLM_Localization.localize(imp, p);
        double[][] d = SMLM_Localization.driftFromLocalizations(locs, W, H, nF, p);
        // correction at the last segment centre should be about (-300, +200) nm x (frame fraction)
        int f = nF - 5; double t = f / (double) (nF - 1);
        System.out.printf(java.util.Locale.ROOT, "drift correction at frame %d: dx %.1f nm (expected %.1f), dy %.1f nm (expected %.1f)%n",
                f, d[f][0] * pix, -300 * t + 300 * 4.5 / 99, d[f][1] * pix, 200 * t - 200 * 4.5 / 99);
        // the first segment (frames 0-9, mean drift at t = 4.5/99) is the reference
        assertEquals(-300 * t + 300 * 4.5 / 99, d[f][0] * pix, 15);
        assertEquals(200 * t - 200 * 4.5 / 99, d[f][1] * pix, 15);
    }

    @Test
    public void mortensenLimit() {
        // no background, no pixelation limit: sigma / sqrt(N) * 4/3 (x sqrt 2 for EMCCD)
        double u = SMLM_Localization.mortensen(10000, 0, 130, 1e-6, false);
        assertEquals(130 / 100.0 * 4 / 3.0, u, 1e-3);
        assertEquals(u * Math.sqrt(2), SMLM_Localization.mortensen(10000, 0, 130, 1e-6, true), 1e-3);
    }
}
