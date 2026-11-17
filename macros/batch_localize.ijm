// batch_localize.ijm - run SMLM Localization on every .tif/.tiff in a folder and save one CSV (+ render) per file.
// Interactive: asks for the folders and the pixel size.
// Headless:    ImageJ-<platform> --headless --console -macro batch_localize.ijm "/input/,/output/,107"
// Edit `options` below for non-default parameters (keys as in the plugin dialog / README).
arg = getArgument();
if (arg != "") {
    a = split(arg, ",");
    inDir = a[0]; outDir = a[1]; pixel = parseFloat(a[2]);
} else {
    inDir = getDirectory("Folder with time series");
    outDir = getDirectory("Output folder");
    pixel = getNumber("Pixel size (nm)", 107);
}
if (!endsWith(inDir, File.separator)) inDir = inDir + File.separator;
if (!endsWith(outDir, File.separator)) outDir = outDir + File.separator;
options = "channel=1 pixel_size=" + pixel + " dog_sigma1=107 dog_sigma2=268 threshold=4 min_separation=321 " +
          "fit_window=960 sigma_min=86 sigma_max=235 camera_offset=0 photons_per_adu=0 emccd " +
          "render render_pixel=10 render_blur=15 time_bins=0";
setBatchMode(true);
list = getFileList(inDir);
n = 0;
for (i = 0; i < list.length; i++) {
    name = list[i];
    if (!(endsWith(toLowerCase(name), ".tif") || endsWith(toLowerCase(name), ".tiff"))) continue;
    open(inDir + name);
    base = File.nameWithoutExtension;
    run("SMLM Localization...", options + " save_csv=[" + outDir + base + "_localisations.csv]");
    close("*");
    n++;
}
print("batch_localize: " + n + " file(s) processed, results in " + outDir);
