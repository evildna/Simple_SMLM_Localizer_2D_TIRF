/*
 * SMLM_Localization - single-molecule localisation for sparse-emitter time series (ImageJ / Fiji plugin).
 * Copyright (c) 2026 Xian Hu / University of Oslo. MIT License, see LICENSE.
 */
package no.uio.imaging;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.GenericDialog;
import ij.measure.Calibration;
import ij.measure.ResultsTable;
import ij.plugin.PlugIn;
import ij.plugin.filter.GaussianBlur;
import ij.process.FHT;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import java.awt.GraphicsEnvironment;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Locale;

/**
 * Single-molecule localisation for sparse-emitter time series (live-cell TIRF, PALM, STORM, PAINT).
 *
 * <p>Per frame of the selected channel:
 * <ol>
 *   <li>difference-of-Gaussians band-pass (default sigma 107 / 268 nm);</li>
 *   <li>candidates: local maxima of the filtered image within a (2m+1)^2 neighbourhood (m = minimum separation, default
 *       321 nm, i.e. 3 px at 107 nm/px and 2 px at 160 nm/px) exceeding k x robust SD (1.4826 x MAD) of the filtered
 *       image (k = 4);</li>
 *   <li>least-squares fit of a symmetric 2D Gaussian plus constant background in a (2r+1)^2 window of the raw frame
 *       (r = round(fit window / pixel / 2), default window 960 nm), bounded Levenberg-Marquardt;</li>
 *   <li>acceptance: amplitude &gt; 0 and sigma inside the gate (default 86-235 nm).</li>
 * </ol>
 * Every spatial parameter is defined in nm and converted with the pixel size, so data at different magnifications
 * are analysed identically. The defaults reproduce the Python analysis of the accompanying manuscript.
 *
 * <p>Optional: residual-drift correction from the localisations (cross-correlation of localisation histograms of
 * successive segments with the first segment), a super-resolution density render and a time-binned render stack.
 *
 * <p>Output: a ResultsTable whose column names follow ThunderSTORM's conventions (File &gt; Save As... writes CSV).
 *
 * <p>The source is written in the subset of Java that Groovy 2.4+ also compiles; the Groovy script
 * SMLM_Localization_Fiji.groovy is generated from this file by build.sh (do not edit it by hand).
 */
public class SMLM_Localization implements PlugIn {

    public static final String VERSION = "1.0.0";
    public static final String URL = "https://github.com/evildna/Simple_SMLM_Localizer_2D_TIRF";

    /** Fit constants of the reference implementation, in nm (= 0.6, 1.3 and 4 x 107 nm); not user parameters. */
    static final double FIT_SIGMA_LOWER_NM = 64.2;
    static final double FIT_SIGMA_INIT_NM = 139.1;
    static final double FIT_SIGMA_UPPER_NM = 428.0;
    /** Bound on the fitted centre relative to the candidate pixel (px), as in the reference implementation. */
    static final double FIT_SHIFT_MAX_PX = 2.0;

    /** All analysis parameters. Spatial parameters are in nm; the pixel equivalents are derived. */
    public static class Params {
        public int channel = 1;
        public double pixelNm = 107.0;
        public double dogSigma1Nm = 107.0;
        public double dogSigma2Nm = 268.0;
        public double thresholdK = 4.0;
        public double minSeparationNm = 321.0;
        public double fitWindowNm = 960.0;
        public double sigmaMinNm = 86.0;
        public double sigmaMaxNm = 235.0;
        public double offsetAdu = 0.0;
        public double photonsPerAdu = 0.0;      // 0 = not calibrated: intensities stay in ADU, no uncertainty column
        public boolean emccd = true;
        public boolean driftCorrection = false;
        public int driftSegment = 25;
        public int driftUpsampling = 2;
        public boolean render = true;
        public double renderPixelNm = 10.0;
        public double renderBlurNm = 15.0;
        public int timeBins = 0;               // 0 = no time-binned stack
        public String saveCsv = "";            // results table path (renders are saved next to it); empty = none
        public int maxIterations = 200;

        public double dogSigma1Px() { return dogSigma1Nm / pixelNm; }
        public double dogSigma2Px() { return dogSigma2Nm / pixelNm; }
        public double sigmaMinPx() { return sigmaMinNm / pixelNm; }
        public double sigmaMaxPx() { return sigmaMaxNm / pixelNm; }
        public int fitRadiusPx() { return Math.max(2, (int) Math.round(fitWindowNm / pixelNm / 2.0)); }
        public int minSeparationPx() { return Math.max(2, (int) Math.round(minSeparationNm / pixelNm)); }
        public boolean calibrated() { return photonsPerAdu > 0; }

        public String describe() {
            int r = fitRadiusPx();
            int m = minSeparationPx();
            return String.format(Locale.ROOT,
                "pixel %.2f nm | DoG sigma %.1f / %.1f nm (%.3f / %.3f px) | threshold %.2f x robust SD | " +
                "neighbourhood %dx%d px (min. separation %.0f nm) | fit window %dx%d px (%.0f nm requested) | " +
                "sigma gate %.1f-%.1f nm (%.3f-%.3f px) | offset %.1f ADU | photons/ADU %s | EMCCD %s",
                pixelNm, dogSigma1Nm, dogSigma2Nm, dogSigma1Px(), dogSigma2Px(), thresholdK,
                2 * m + 1, 2 * m + 1, minSeparationNm, 2 * r + 1, 2 * r + 1, fitWindowNm,
                sigmaMinNm, sigmaMaxNm, sigmaMinPx(), sigmaMaxPx(), offsetAdu,
                calibrated() ? String.format(Locale.ROOT, "%.5f", photonsPerAdu) : "not set (ADU)", emccd ? "yes" : "no");
        }
    }

    /** One localisation. x, y in camera pixels with the centre of the first pixel at 0; sigma in px. */
    public static class Loc {
        public int frame;          // 0-based
        public double x, y, sigma, amplitude, background, chi2;
        public double intensity, offset, uncertainty;   // ADU or photons; uncertainty in nm (NaN if not calibrated)
    }

    // ================================================================== ImageJ entry point
    public void run(String arg) {
        ImagePlus imp = IJ.getImage();
        if (imp == null) return;
        Params p = new Params();
        Calibration cal = imp.getCalibration();
        double calNm = calibrationNm(cal);
        if (calNm > 0) p.pixelNm = calNm;
        if (!showDialog(imp, p)) return;

        IJ.log("SMLM_Localization " + VERSION + " on '" + imp.getTitle() + "', channel " + p.channel);
        IJ.log("  " + p.describe());
        long t0 = System.currentTimeMillis();
        List<Loc> locs = localize(imp, p);
        int nF = frameCount(imp);
        IJ.log(String.format(Locale.ROOT, "  %d localisations in %d frames (%.1f per frame), %.1f s",
                locs.size(), nF, locs.size() / (double) Math.max(1, nF), (System.currentTimeMillis() - t0) / 1000.0));
        if (!p.calibrated())
            IJ.log("  photons per ADU not set: intensity and offset are in ADU; no localisation uncertainty is reported.");

        double[][] drift = null;
        if (p.driftCorrection) {
            if (nF >= 2 * p.driftSegment) {
                drift = driftFromLocalizations(locs, imp.getWidth(), imp.getHeight(), nF, p);
                applyDrift(locs, drift);
                IJ.log(String.format(Locale.ROOT, "  localisation-based drift at last frame: dx = %.3f, dy = %.3f px (%.1f, %.1f nm)",
                        drift[nF - 1][0], drift[nF - 1][1], drift[nF - 1][0] * p.pixelNm, drift[nF - 1][1] * p.pixelNm));
            } else {
                IJ.log("  drift correction skipped: fewer than 2 segments of " + p.driftSegment + " frames.");
            }
        }

        ResultsTable rt = toResultsTable(locs, drift, p, frameIntervalSeconds(cal));
        String base = p.saveCsv.length() == 0 ? null : (p.saveCsv.toLowerCase(Locale.ROOT).endsWith(".csv")
                ? p.saveCsv.substring(0, p.saveCsv.length() - 4) : p.saveCsv);
        if (base != null) {
            rt.save(base + ".csv");
            IJ.log("  saved " + base + ".csv");
        }
        boolean gui = !GraphicsEnvironment.isHeadless();   // older ImageJ throws on show() when headless
        if (gui) rt.show("Localizations_" + imp.getTitle());
        if (p.render) {
            ImagePlus r = render(locs, imp.getWidth(), imp.getHeight(), p, imp.getTitle());
            if (base != null) IJ.saveAsTiff(r, base + "_SMLM_" + fmt(p.renderPixelNm) + "nm.tif");
            if (gui) r.show();
        }
        if (p.timeBins > 0) {
            ImagePlus tb = renderTimeBins(locs, imp.getWidth(), imp.getHeight(), nF, p, imp.getTitle());
            if (base != null) IJ.saveAsTiff(tb, base + "_SMLM_timebins_" + fmt(p.renderPixelNm) + "nm.tif");
            if (gui) tb.show();
        }
    }

    /** Parameter dialog; macro keys are the first word of each label (e.g. dog_sigma1=107). */
    public static boolean showDialog(ImagePlus imp, Params p) {
        int nC = imp.getNChannels();
        GenericDialog gd = new GenericDialog("SMLM Localization " + VERSION);
        gd.addMessage("Spatial parameters in nm (converted with the pixel size)");
        gd.addNumericField("Channel (1-" + nC + ")", p.channel, 0);
        gd.addNumericField("Pixel_size (nm)", p.pixelNm, 2);
        gd.addNumericField("DoG_sigma1 (nm)", p.dogSigma1Nm, 1);
        gd.addNumericField("DoG_sigma2 (nm)", p.dogSigma2Nm, 1);
        gd.addNumericField("Threshold (k x robust SD)", p.thresholdK, 2);
        gd.addNumericField("Min_separation (nm)", p.minSeparationNm, 1);
        gd.addNumericField("Fit_window (nm)", p.fitWindowNm, 1);
        gd.addNumericField("Sigma_min (nm)", p.sigmaMinNm, 1);
        gd.addNumericField("Sigma_max (nm)", p.sigmaMaxNm, 1);
        gd.addMessage("Camera (photons per ADU = e-/ADU / EM gain; 0 = keep ADU)");
        gd.addNumericField("Camera_offset (ADU)", p.offsetAdu, 1);
        gd.addNumericField("Photons_per_ADU", p.photonsPerAdu, 5);
        gd.addCheckbox("EMCCD (excess-noise factor 2 in uncertainty)", p.emccd);
        gd.addMessage("Optional");
        gd.addCheckbox("Drift_correction from localisations", p.driftCorrection);
        gd.addNumericField("Drift_segment (frames)", p.driftSegment, 0);
        gd.addNumericField("Drift_upsampling (x camera pixel)", p.driftUpsampling, 0);
        gd.addCheckbox("Render super-resolution image", p.render);
        gd.addNumericField("Render_pixel (nm)", p.renderPixelNm, 1);
        gd.addNumericField("Render_blur (nm)", p.renderBlurNm, 1);
        gd.addNumericField("Time_bins (0 = none)", p.timeBins, 0);
        gd.addStringField("Save_CSV (path; empty = do not save)", p.saveCsv, 30);
        gd.addHelp(URL);
        gd.showDialog();
        if (gd.wasCanceled()) return false;
        p.channel = Math.min(nC, Math.max(1, (int) gd.getNextNumber()));
        p.pixelNm = gd.getNextNumber();
        p.dogSigma1Nm = gd.getNextNumber();
        p.dogSigma2Nm = gd.getNextNumber();
        p.thresholdK = gd.getNextNumber();
        p.minSeparationNm = gd.getNextNumber();
        p.fitWindowNm = gd.getNextNumber();
        p.sigmaMinNm = gd.getNextNumber();
        p.sigmaMaxNm = gd.getNextNumber();
        p.offsetAdu = gd.getNextNumber();
        p.photonsPerAdu = Math.max(0, gd.getNextNumber());
        p.emccd = gd.getNextBoolean();
        p.driftCorrection = gd.getNextBoolean();
        p.driftSegment = Math.max(2, (int) gd.getNextNumber());
        p.driftUpsampling = Math.max(1, (int) gd.getNextNumber());
        p.render = gd.getNextBoolean();
        p.renderPixelNm = Math.max(1, gd.getNextNumber());
        p.renderBlurNm = Math.max(0, gd.getNextNumber());
        p.timeBins = Math.max(0, (int) gd.getNextNumber());
        p.saveCsv = gd.getNextString().trim();
        if (!(p.pixelNm > 0) || !(p.dogSigma2Nm > p.dogSigma1Nm) || !(p.sigmaMaxNm > p.sigmaMinNm)) {
            IJ.error("SMLM Localization", "Need pixel size > 0, DoG sigma2 > sigma1 and sigma max > sigma min.");
            return false;
        }
        return true;
    }

    // ================================================================== stack handling
    /** Number of time points: frames of a hyperstack, or slices of a plain stack. */
    public static int frameCount(ImagePlus imp) {
        int nZ = imp.getNSlices();
        int nT = imp.getNFrames();
        if (nT == 1 && nZ > 1) return nZ;
        return Math.max(1, nT);
    }

    /** Pixel size in nm from the image calibration, or -1 if the unit is not a length unit ImageJ reports. */
    public static double calibrationNm(Calibration cal) {
        if (cal == null || !cal.scaled()) return -1;
        String u = cal.getUnit() == null ? "" : cal.getUnit().toLowerCase(Locale.ROOT);
        if (u.equals("nm") || u.startsWith("nano")) return cal.pixelWidth;
        if (u.equals("um") || u.equals("\u00b5m") || u.equals("micron") || u.startsWith("micro")) return cal.pixelWidth * 1000.0;
        return -1;
    }

    /** Frame interval in seconds from the calibration (0 if unknown). */
    public static double frameIntervalSeconds(Calibration cal) {
        if (cal == null || !(cal.frameInterval > 0)) return 0;
        String u = cal.getTimeUnit() == null ? "sec" : cal.getTimeUnit().toLowerCase(Locale.ROOT);
        if (u.startsWith("ms") || u.startsWith("milli")) return cal.frameInterval / 1000.0;
        if (u.startsWith("min")) return cal.frameInterval * 60.0;
        if (u.startsWith("h")) return cal.frameInterval * 3600.0;
        return cal.frameInterval;
    }

    /** Localise all time points of the selected channel (central slice of z-stacks). */
    public static List<Loc> localize(ImagePlus imp, Params p) {
        List<Loc> locs = new ArrayList<Loc>();
        ImageStack st = imp.getStack();
        int nZ = imp.getNSlices();
        int nT = imp.getNFrames();
        boolean framesAreSlices = (nT == 1 && nZ > 1);
        int nF = frameCount(imp);
        int c = Math.min(imp.getNChannels(), Math.max(1, p.channel));
        for (int f = 0; f < nF; f++) {
            if (f % 10 == 0) {
                IJ.showStatus("SMLM: frame " + (f + 1) + "/" + nF + " (" + locs.size() + " localisations)");
                IJ.showProgress(f, nF);
            }
            int idx = framesAreSlices ? imp.getStackIndex(c, f + 1, 1) : imp.getStackIndex(c, (int) ((nZ + 1) / 2), f + 1);
            FloatProcessor fp = st.getProcessor(idx).convertToFloatProcessor();
            localizeFrame((float[]) fp.getPixels(), fp.getWidth(), fp.getHeight(), f, p, locs);
        }
        IJ.showProgress(1.0);
        return locs;
    }

    // ================================================================== detection and fitting
    /** Detect and fit one frame (raw pixel values, row-major) and append accepted localisations to out. */
    public static void localizeFrame(float[] raw, int W, int H, int frame, Params p, List<Loc> out) {
        float[] g1 = gaussianFilter(raw, W, H, p.dogSigma1Px());
        float[] g2 = gaussianFilter(raw, W, H, p.dogSigma2Px());
        int n = W * H;
        float[] dog = new float[n];
        for (int i = 0; i < n; i++) dog[i] = g1[i] - g2[i];
        double robustSd = 1.4826 * mad(dog);
        double thr = p.thresholdK * robustSd;
        int r = p.fitRadiusPx();
        int m = p.minSeparationPx();
        int border = r + 1;
        double sLo = FIT_SIGMA_LOWER_NM / p.pixelNm;
        double sHi = FIT_SIGMA_UPPER_NM / p.pixelNm;
        double s0 = FIT_SIGMA_INIT_NM / p.pixelNm;
        double gateLo = p.sigmaMinPx();
        double gateHi = p.sigmaMaxPx();
        for (int y = border; y < H - border; y++) {
            for (int x = border; x < W - border; x++) {
                float v = dog[y * W + x];
                if (!(v > thr)) continue;
                if (!isLocalMax(dog, W, H, x, y, m)) continue;
                double[] q = fitGaussian(raw, W, x, y, r, s0, sLo, sHi, p.maxIterations);
                double amp = q[0];
                double sg = q[3];
                if (!(amp > 0) || !(sg > gateLo) || !(sg < gateHi)) continue;
                Loc l = new Loc();
                l.frame = frame;
                l.x = x + q[1];
                l.y = y + q[2];
                l.sigma = sg;
                l.amplitude = amp;
                l.background = q[4];
                l.chi2 = q[5];
                double nAdu = 2.0 * Math.PI * amp * sg * sg;
                if (p.calibrated()) {
                    l.intensity = nAdu * p.photonsPerAdu;
                    l.offset = Math.max(0, q[4] - p.offsetAdu) * p.photonsPerAdu;
                    l.uncertainty = mortensen(l.intensity, l.offset, sg * p.pixelNm, p.pixelNm, p.emccd);
                } else {
                    l.intensity = nAdu;
                    l.offset = q[4];
                    l.uncertainty = Double.NaN;
                }
                out.add(l);
            }
        }
    }

    /** Maximum within the (2m+1)^2 neighbourhood; ties are resolved in favour of the first pixel in raster order. */
    static boolean isLocalMax(float[] img, int W, int H, int x, int y, int m) {
        float v = img[y * W + x];
        int y0 = Math.max(0, y - m);
        int y1 = Math.min(H - 1, y + m);
        int x0 = Math.max(0, x - m);
        int x1 = Math.min(W - 1, x + m);
        for (int yy = y0; yy <= y1; yy++) {
            for (int xx = x0; xx <= x1; xx++) {
                if (xx == x && yy == y) continue;
                float u = img[yy * W + xx];
                if (u > v) return false;
                if (u == v && (yy < y || (yy == y && xx < x))) return false;
            }
        }
        return true;
    }

    /**
     * Localisation uncertainty (nm), Mortensen et al. 2010 least-squares form with pixelation term; N and b in photons
     * (b = background photons per pixel), s and a (pixel) in nm. EMCCD doubles the variance (excess-noise factor 2).
     */
    public static double mortensen(double N, double b, double sNm, double aNm, boolean emccd) {
        if (!(N > 0)) return Double.NaN;
        double sa2 = sNm * sNm + aNm * aNm / 12.0;
        double var = sa2 / N * (16.0 / 9.0 + 8.0 * Math.PI * sa2 * b / (N * aNm * aNm));
        if (emccd) var *= 2.0;
        return Math.sqrt(var);
    }

    /** Median absolute deviation from the median (numpy convention for even counts). */
    static double mad(float[] a) {
        float[] s = Arrays.copyOf(a, a.length);
        Arrays.sort(s);
        double med = medianSorted(s);
        for (int i = 0; i < s.length; i++) s[i] = (float) Math.abs(s[i] - med);
        Arrays.sort(s);
        return medianSorted(s);
    }

    static double medianSorted(float[] s) {
        int n = s.length;
        if (n == 0) return Double.NaN;
        int h = (int) (n / 2);
        if (n % 2 == 1) return s[h];
        return 0.5 * ((double) s[h - 1] + (double) s[h]);
    }

    /**
     * Separable Gaussian filter, kernel truncated at 4 sigma, mirror ("reflect", d c b a | a b c d | d c b a)
     * boundaries; the same definition as scipy.ndimage.gaussian_filter.
     */
    public static float[] gaussianFilter(float[] img, int W, int H, double sigma) {
        int rad = (int) (4.0 * sigma + 0.5);
        double[] k = new double[2 * rad + 1];
        double sum = 0;
        for (int i = -rad; i <= rad; i++) {
            k[i + rad] = Math.exp(-0.5 * i * i / (sigma * sigma));
            sum += k[i + rad];
        }
        for (int i = 0; i < k.length; i++) k[i] /= sum;
        float[] tmp = new float[W * H];
        float[] out = new float[W * H];
        for (int y = 0; y < H; y++) {
            int row = y * W;
            for (int x = 0; x < W; x++) {
                double acc = 0;
                for (int i = -rad; i <= rad; i++) acc += k[i + rad] * img[row + reflect(x + i, W)];
                tmp[row + x] = (float) acc;
            }
        }
        for (int y = 0; y < H; y++) {
            for (int x = 0; x < W; x++) {
                double acc = 0;
                for (int i = -rad; i <= rad; i++) acc += k[i + rad] * tmp[reflect(y + i, H) * W + x];
                out[y * W + x] = (float) acc;
            }
        }
        return out;
    }

    static int reflect(int i, int n) {
        if (n == 1) return 0;
        int period = 2 * n;
        int j = i % period;
        if (j < 0) j += period;
        return j < n ? j : period - 1 - j;
    }

    /**
     * Bounded Levenberg-Marquardt fit of z = A exp(-((x-x0)^2 + (y-y0)^2) / (2 s^2)) + b in the (2r+1)^2 window centred on
     * (cx, cy). Start values and bounds as in the reference implementation: b0 = 20th percentile of the window,
     * A0 = max - b0, s0 given; A, b &gt;= 0, |x0|, |y0| &lt;= 2 px, s in [sLo, sHi].
     * Returns {A, x0, y0, s, b, reduced chi2}.
     */
    static double[] fitGaussian(float[] px, int W, int cx, int cy, int r, double s0, double sLo, double sHi, int maxIter) {
        int side = 2 * r + 1;
        int n = side * side;
        double[] z = new double[n];
        int k = 0;
        double mx = -Double.MAX_VALUE;
        for (int dy = -r; dy <= r; dy++) {
            for (int dx = -r; dx <= r; dx++) {
                double v = px[(cy + dy) * W + cx + dx];
                z[k++] = v;
                if (v > mx) mx = v;
            }
        }
        double[] zs = Arrays.copyOf(z, n);
        Arrays.sort(zs);
        double pos = 0.2 * (n - 1);
        int i0 = (int) Math.floor(pos);
        double b0 = zs[i0] + (pos - i0) * (zs[Math.min(n - 1, i0 + 1)] - zs[i0]);
        double[] lo = new double[5];
        double[] hi = new double[5];
        lo[0] = 0; hi[0] = Double.POSITIVE_INFINITY;
        lo[1] = -FIT_SHIFT_MAX_PX; hi[1] = FIT_SHIFT_MAX_PX;
        lo[2] = -FIT_SHIFT_MAX_PX; hi[2] = FIT_SHIFT_MAX_PX;
        lo[3] = sLo; hi[3] = sHi;
        lo[4] = 0; hi[4] = Double.POSITIVE_INFINITY;
        double[] p = new double[5];
        p[0] = mx - b0; p[1] = 0; p[2] = 0; p[3] = Math.min(sHi, Math.max(sLo, s0)); p[4] = Math.max(0, b0);
        double lambda = 1e-3;
        double chi = chi2(p, z, r);
        double[] J = new double[5];
        for (int iter = 0; iter < maxIter; iter++) {
            double[][] jtj = new double[5][5];
            double[] jtr = new double[5];
            double s2 = p[3] * p[3];
            k = 0;
            for (int dy = -r; dy <= r; dy++) {
                for (int dx = -r; dx <= r; dx++) {
                    double ex = dx - p[1];
                    double ey = dy - p[2];
                    double e = Math.exp(-(ex * ex + ey * ey) / (2 * s2));
                    double res = z[k++] - (p[0] * e + p[4]);
                    J[0] = e;
                    J[1] = p[0] * e * ex / s2;
                    J[2] = p[0] * e * ey / s2;
                    J[3] = p[0] * e * (ex * ex + ey * ey) / (s2 * p[3]);
                    J[4] = 1.0;
                    for (int a = 0; a < 5; a++) {
                        jtr[a] += J[a] * res;
                        for (int c = 0; c < 5; c++) jtj[a][c] += J[a] * J[c];
                    }
                }
            }
            // active set: parameters at a bound whose descent direction points outward are held fixed
            boolean[] fixed = new boolean[5];
            for (int a = 0; a < 5; a++) {
                if (jtj[a][a] <= 1e-30) fixed[a] = true;
                if (p[a] <= lo[a] && jtr[a] < 0) fixed[a] = true;
                if (p[a] >= hi[a] && jtr[a] > 0) fixed[a] = true;
            }
            double[][] mat = new double[5][5];
            double[] rhs = new double[5];
            for (int a = 0; a < 5; a++) {
                for (int c = 0; c < 5; c++) mat[a][c] = (fixed[a] || fixed[c]) ? 0.0d : jtj[a][c];
                mat[a][a] = fixed[a] ? 1.0d : jtj[a][a] * (1.0d + lambda);
                rhs[a] = fixed[a] ? 0.0d : jtr[a];
            }
            double[] d = solve(mat, rhs);
            if (d == null) break;
            double[] q = new double[5];
            double stepMax = 0;
            for (int a = 0; a < 5; a++) {
                q[a] = Math.min(hi[a], Math.max(lo[a], p[a] + d[a]));
                stepMax = Math.max(stepMax, Math.abs(q[a] - p[a]) / (Math.abs(p[a]) + 1e-8));
            }
            double chiq = chi2(q, z, r);
            if (chiq <= chi) {
                boolean converged = (chi - chiq) <= 1e-10 * chi || stepMax < 1e-10;
                p = q;
                chi = chiq;
                lambda = Math.max(lambda * 0.3d, 1e-9d);
                if (converged) break;
            } else {
                lambda *= 10;
                if (lambda > 1e8) break;
            }
        }
        double[] res = new double[6];
        for (int a = 0; a < 5; a++) res[a] = p[a];
        res[5] = chi / (n - 5);
        return res;
    }

    static double chi2(double[] p, double[] z, int r) {
        double s = 0;
        int k = 0;
        double s2 = p[3] * p[3];
        for (int dy = -r; dy <= r; dy++) {
            for (int dx = -r; dx <= r; dx++) {
                double ex = dx - p[1];
                double ey = dy - p[2];
                double res = z[k++] - (p[0] * Math.exp(-(ex * ex + ey * ey) / (2 * s2)) + p[4]);
                s += res * res;
            }
        }
        return s;
    }

    /** Gaussian elimination with partial pivoting; null if singular. */
    static double[] solve(double[][] A, double[] b) {
        int n = b.length;
        double[][] M = new double[n][n + 1];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) M[i][j] = A[i][j];
            M[i][n] = b[i];
        }
        for (int c = 0; c < n; c++) {
            int piv = c;
            for (int i = c + 1; i < n; i++) if (Math.abs(M[i][c]) > Math.abs(M[piv][c])) piv = i;
            if (Math.abs(M[piv][c]) < 1e-300) return null;
            double[] tmp = M[c];
            M[c] = M[piv];
            M[piv] = tmp;
            for (int i = c + 1; i < n; i++) {
                double f = M[i][c] / M[c][c];
                for (int j = c; j <= n; j++) M[i][j] -= f * M[c][j];
            }
        }
        double[] x = new double[n];
        for (int i = n - 1; i >= 0; i--) {
            double s = M[i][n];
            for (int j = i + 1; j < n; j++) s -= M[i][j] * x[j];
            x[i] = s / M[i][i];
        }
        return x;
    }

    // ================================================================== drift from localisations
    /**
     * Residual drift from the localisations: histograms of successive segments (driftSegment frames, driftUpsampling x
     * camera resolution, Gaussian sigma 1 bin) are cross-correlated with the first segment (FHT, 3-point parabolic peak);
     * the sign is verified by re-correlating the most-shifted segment after applying +/- the shift. Shifts are assigned
     * to the segment centres and interpolated linearly. Returns, per frame, the correction {dx, dy} in pixels that is
     * ADDED to x, y.
     */
    public static double[][] driftFromLocalizations(List<Loc> locs, int W, int H, int nF, Params p) {
        int seg = p.driftSegment;
        int up = p.driftUpsampling;
        int nSeg = (int) (nF / seg);
        int N = pow2(Math.max(W, H) * up);
        FloatProcessor[] hist = new FloatProcessor[nSeg];
        for (int s = 0; s < nSeg; s++) hist[s] = new FloatProcessor(N, N);
        for (int i = 0; i < locs.size(); i++) {
            Loc l = locs.get(i);
            int s = (int) (l.frame / seg);
            if (s >= nSeg) continue;
            int xr = (int) Math.floor((l.x + 0.5) * up);
            int yr = (int) Math.floor((l.y + 0.5) * up);
            if (xr >= 0 && yr >= 0 && xr < N && yr < N) ((float[]) hist[s].getPixels())[yr * N + xr] += 1f;
        }
        GaussianBlur gb = new GaussianBlur();
        for (int s = 0; s < nSeg; s++) gb.blurGaussian(hist[s], 1.0, 1.0, 0.002);
        FHT ref = new FHT(hist[0]);
        ref.transform();
        double[] sx = new double[nSeg];
        double[] sy = new double[nSeg];
        int sMax = 0;
        for (int s = 1; s < nSeg; s++) {
            IJ.showStatus("SMLM drift: segment " + s + "/" + (nSeg - 1));
            double[] sh = peakShift(ref, hist[s]);
            sx[s] = sh[0] / up;
            sy[s] = sh[1] / up;
            if (Math.hypot(sx[s], sy[s]) > Math.hypot(sx[sMax], sy[sMax])) sMax = s;
        }
        if (nSeg > 1 && Math.hypot(sx[sMax], sy[sMax]) * up > 0.5) {
            double rPlus = residualAfter(ref, hist[sMax], sx[sMax] * up, sy[sMax] * up);
            double rMinus = residualAfter(ref, hist[sMax], -sx[sMax] * up, -sy[sMax] * up);
            if (rMinus < rPlus) {
                for (int s = 0; s < nSeg; s++) { sx[s] = -sx[s]; sy[s] = -sy[s]; }
                IJ.log("  drift: shift sign inverted by self-test");
            }
        }
        double[][] out = new double[nF][2];
        for (int f = 0; f < nF; f++) {
            double tb = (f - (seg - 1) / 2.0) / seg;
            int b0 = (int) Math.floor(tb);
            double w = tb - b0;
            if (b0 < 0) { out[f][0] = sx[0]; out[f][1] = sy[0]; }
            else if (b0 >= nSeg - 1) { out[f][0] = sx[nSeg - 1]; out[f][1] = sy[nSeg - 1]; }
            else {
                out[f][0] = sx[b0] * (1 - w) + sx[b0 + 1] * w;
                out[f][1] = sy[b0] * (1 - w) + sy[b0 + 1] * w;
            }
        }
        return out;
    }

    public static void applyDrift(List<Loc> locs, double[][] drift) {
        for (int i = 0; i < locs.size(); i++) {
            Loc l = locs.get(i);
            l.x += drift[l.frame][0];
            l.y += drift[l.frame][1];
        }
    }

    static int pow2(int n) {
        int q = 1;
        while (q < n) q = q * 2;
        return q;
    }

    /** Position of the cross-correlation peak of ref (transformed) with mov, relative to the centre (pixels). */
    static double[] peakShift(FHT ref, FloatProcessor mov) {
        FHT h = new FHT(mov);
        h.transform();
        FHT cc = ref.conjugateMultiply(h);
        cc.inverseTransform();
        cc.swapQuadrants();
        int N = cc.getWidth();
        float[] c = (float[]) cc.getPixels();
        int im = 0;
        for (int i = 1; i < c.length; i++) if (c[i] > c[im]) im = i;
        int py = (int) (im / N);
        int px = im - py * N;
        double sx = px;
        double sy = py;
        if (px > 0 && px < N - 1) {
            double d = c[im - 1] - 2 * c[im] + c[im + 1];
            if (d != 0) sx += 0.5 * (c[im - 1] - c[im + 1]) / d;
        }
        if (py > 0 && py < N - 1) {
            double d = c[im - N] - 2 * c[im] + c[im + N];
            if (d != 0) sy += 0.5 * (c[im - N] - c[im + N]) / d;
        }
        int half = (int) (N / 2);
        double[] out = new double[2];
        out[0] = sx - half;
        out[1] = sy - half;
        return out;
    }

    static double residualAfter(FHT ref, FloatProcessor mov, double sx, double sy) {
        FloatProcessor t = (FloatProcessor) mov.duplicate();
        t.setInterpolationMethod(ImageProcessor.BILINEAR);
        t.setBackgroundValue(0);
        t.translate(sx, sy);
        double[] r = peakShift(ref, t);
        return Math.hypot(r[0], r[1]);
    }

    // ================================================================== output
    /**
     * Results table with ThunderSTORM-style column names. Coordinates use ThunderSTORM's convention (the top-left corner
     * of the image is 0 nm, so the centre of the first pixel is at 0.5 px).
     */
    public static ResultsTable toResultsTable(List<Loc> locs, double[][] drift, Params p, double frameIntervalS) {
        ResultsTable rt = new ResultsTable();
        rt.setPrecision(4);
        String iu = p.calibrated() ? "photon" : "ADU";
        for (int i = 0; i < locs.size(); i++) {
            Loc l = locs.get(i);
            rt.incrementCounter();
            rt.addValue("id", i + 1);
            rt.addValue("frame", l.frame + 1);
            if (frameIntervalS > 0) rt.addValue("t [s]", l.frame * frameIntervalS);
            rt.addValue("x [nm]", (l.x + 0.5) * p.pixelNm);
            rt.addValue("y [nm]", (l.y + 0.5) * p.pixelNm);
            rt.addValue("sigma [nm]", l.sigma * p.pixelNm);
            rt.addValue("intensity [" + iu + "]", l.intensity);
            rt.addValue("offset [" + iu + "]", l.offset);
            if (p.calibrated()) rt.addValue("uncertainty [nm]", l.uncertainty);
            rt.addValue("chi2", l.chi2);
            if (drift != null) {
                rt.addValue("drift_x [nm]", drift[l.frame][0] * p.pixelNm);
                rt.addValue("drift_y [nm]", drift[l.frame][1] * p.pixelNm);
            }
        }
        return rt;
    }

    /** Density render: localisation histogram at renderPixelNm, Gaussian blur renderBlurNm; calibrated in nm. */
    public static ImagePlus render(List<Loc> locs, int W, int H, Params p, String title) {
        int RW = (int) Math.ceil(W * p.pixelNm / p.renderPixelNm);
        int RH = (int) Math.ceil(H * p.pixelNm / p.renderPixelNm);
        float[] img = histogram(locs, RW, RH, p, 0, Integer.MAX_VALUE);
        if (p.renderBlurNm > 0) img = gaussianFilter(img, RW, RH, p.renderBlurNm / p.renderPixelNm);
        ImagePlus out = new ImagePlus(title + "_SMLM_" + fmt(p.renderPixelNm) + "nm", new FloatProcessor(RW, RH, img, null));
        calibrateNm(out, p.renderPixelNm);
        out.getProcessor().resetMinAndMax();
        return out;
    }

    /** Stack of nBins renders of equal frame ranges (colour-code with Image > Hyperstacks > Temporal-Color Code). */
    public static ImagePlus renderTimeBins(List<Loc> locs, int W, int H, int nF, Params p, String title) {
        int RW = (int) Math.ceil(W * p.pixelNm / p.renderPixelNm);
        int RH = (int) Math.ceil(H * p.pixelNm / p.renderPixelNm);
        ImageStack st = new ImageStack(RW, RH);
        int nb = Math.max(1, p.timeBins);
        for (int b = 0; b < nb; b++) {
            int f0 = (int) Math.floor(b * (double) nF / nb);
            int f1 = (int) Math.floor((b + 1) * (double) nF / nb);
            float[] img = histogram(locs, RW, RH, p, f0, f1);
            if (p.renderBlurNm > 0) img = gaussianFilter(img, RW, RH, p.renderBlurNm / p.renderPixelNm);
            st.addSlice("frames " + (f0 + 1) + "-" + f1, new FloatProcessor(RW, RH, img, null));
        }
        ImagePlus out = new ImagePlus(title + "_SMLM_timebins_" + fmt(p.renderPixelNm) + "nm", st);
        out.setDimensions(1, 1, nb);
        out.setOpenAsHyperStack(true);
        calibrateNm(out, p.renderPixelNm);
        return out;
    }

    static float[] histogram(List<Loc> locs, int RW, int RH, Params p, int f0, int f1) {
        float[] img = new float[RW * RH];
        double sc = p.pixelNm / p.renderPixelNm;
        for (int i = 0; i < locs.size(); i++) {
            Loc l = locs.get(i);
            if (l.frame < f0 || l.frame >= f1) continue;
            int xr = (int) Math.floor((l.x + 0.5) * sc);
            int yr = (int) Math.floor((l.y + 0.5) * sc);
            if (xr >= 0 && yr >= 0 && xr < RW && yr < RH) img[yr * RW + xr] += 1f;
        }
        return img;
    }

    static void calibrateNm(ImagePlus imp, double nm) {
        Calibration c = imp.getCalibration();
        c.pixelWidth = nm;
        c.pixelHeight = nm;
        c.setUnit("nm");
    }

    static String fmt(double v) {
        return (v == Math.rint(v)) ? String.valueOf((long) v) : String.format(Locale.ROOT, "%.1f", v);
    }
}
