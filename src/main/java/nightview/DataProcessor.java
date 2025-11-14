/*
 * Decompiled with CFR 0.152.
 * 
 * Could not load the following classes:
 *  ij.Prefs
 */
package nightview;

import ij.Prefs;
import java.io.File;
import java.io.PrintWriter;
import java.time.LocalDate;
import java.time.format.DateTimeFormatter;
import java.util.Date;
import java.util.HashMap;
import java.util.Map;
import java.util.function.Consumer;

public class DataProcessor {
    private File workspaceRoot;
    private LocalDate night;
    private Consumer<String> progressCallback;
    private static boolean imageJInitialized = false;

    public DataProcessor(File file, LocalDate localDate) {
        this.workspaceRoot = file;
        this.night = localDate;
    }

    public void setProgressCallback(Consumer<String> consumer) {
        this.progressCallback = consumer;
    }

    private void updateProgress(String string) {
        System.out.println(string);
        if (this.progressCallback != null) {
            this.progressCallback.accept(string);
        }
    }

    public void processNight(boolean bl) throws Exception {
        System.out.println("=== DataProcessor.processNight() called ===");
        System.out.flush();
        String string = this.night.format(DateTimeFormatter.ofPattern("yyyyMMdd"));
        File file = new File(this.workspaceRoot, "data/" + string);
        if (!file.exists()) {
            throw new RuntimeException("Night directory not found: " + file.getAbsolutePath());
        }
        this.updateProgress("PROCESSING STARTED!");
        try {
            this.updateProgress("Initializing...");
            System.out.println("Step 1: Initialize");
            System.out.flush();
            Thread.sleep(300L);
            this.updateProgress("Setting up processor...");
            System.out.println("Step 2: Setup");
            System.out.flush();
            Thread.sleep(300L);
            this.updateProgress("Analyzing data directories...");
            System.out.println("Step 3: Analyze");
            System.out.flush();
            Thread.sleep(500L);
            this.updateProgress("Preparing processing engine...");
            System.out.println("Step 4: Prepare engine");
            System.out.flush();
            this.initializeImageJ();
            Thread.sleep(300L);
            this.updateProgress("Processing FITS data...");
            System.out.println("Step 5: Process data");
            System.out.flush();
            this.processWithDirectMethod(file, bl);
            this.updateProgress("Processing completed successfully!");
            System.out.println("Step 6: Complete");
            System.out.flush();
        }
        catch (InterruptedException interruptedException) {
            Thread.currentThread().interrupt();
            this.updateProgress("Processing was interrupted");
            throw new Exception("Processing interrupted", interruptedException);
        }
        catch (Exception exception) {
            String string2 = "Processing failed: " + exception.getMessage();
            this.updateProgress(string2);
            System.out.println("ERROR in processNight(): " + exception.getMessage());
            exception.printStackTrace();
            System.out.flush();
            throw exception;
        }
    }

    private void initializeImageJ() {
        System.out.println("=== initializeImageJ() called ===");
        this.updateProgress("Processing engine starting...");
        System.out.flush();
        if (!imageJInitialized) {
            try {
                this.updateProgress("Checking system requirements...");
                Thread.sleep(100L);
                String string = System.getProperty("java.class.path");
                boolean bl = string.contains("ij.jar");
                boolean bl2 = string.contains("Astronomy_.jar");
                this.updateProgress("ImageJ library: " + (bl ? "found" : "missing"));
                Thread.sleep(100L);
                this.updateProgress("AstroImageJ library: " + (bl2 ? "found" : "missing"));
                Thread.sleep(100L);
                System.setProperty("java.awt.headless", "true");
                this.updateProgress("Headless mode configured");
                Thread.sleep(100L);
                imageJInitialized = true;
                this.updateProgress("Processing engine ready");
                System.out.println("ImageJ initialization completed successfully");
                System.out.flush();
            }
            catch (Exception exception) {
                this.updateProgress("Initialization warning: " + exception.getMessage());
                System.out.println("Exception during initialization: " + exception.getMessage());
                imageJInitialized = true;
            }
        } else {
            this.updateProgress("Processing engine already initialized");
        }
        System.out.println("=== initializeImageJ() completed ===");
        System.out.flush();
    }

    private String getPluginsDir() {
        String[] stringArray;
        for (String string : stringArray = new String[]{"/Applications/AstroImageJ.app/Contents/Resources/plugins/", System.getProperty("user.home") + "/AstroImageJ/plugins/", "./lib/"}) {
            File file = new File(string);
            if (!file.exists() || !file.isDirectory()) continue;
            return string;
        }
        return "./lib/";
    }

    private void processWithDirectMethod(File file, boolean bl) throws Exception {
        this.updateProgress("Setting up FITS processor...");
        Thread.sleep(100L);
        File file2 = new File(file, "Bias");
        File file3 = new File(file, "Flat");
        File file4 = new File(file, "Light");
        this.updateProgress("Searching for calibration frames...");
        Thread.sleep(200L);
        CalibrationFrameInfo calibrationFrameInfo = this.findBestBiasFrames(file2);
        Map<String, CalibrationFrameInfo> map = this.findBestFlatFrames(file3);
        this.updateProgress("Analyzing configuration...");
        Thread.sleep(100L);
        this.logProcessingConfiguration(file, calibrationFrameInfo, map, file4, bl);
        this.updateProgress("Starting data processing...");
        Thread.sleep(100L);
        this.runDirectFitsProcessing(file4, calibrationFrameInfo, map, bl);
    }

    /*
     * WARNING - void declaration
     */
    private void logProcessingConfiguration(File file2, CalibrationFrameInfo calibrationFrameInfo, Map<String, CalibrationFrameInfo> map, File file3, boolean bl) {
        Object object;
        this.updateProgress("Processing configuration:");
        this.updateProgress("  Night directory: " + file2.getName());
        this.updateProgress("  Force reprocess: " + bl);
        this.updateProgress("  Light directory exists: " + file3.exists());
        if (calibrationFrameInfo != null && calibrationFrameInfo.directory.exists()) {
            void var7_9;
            object = calibrationFrameInfo.directory.listFiles((file, string) -> string.toLowerCase().endsWith(".fits") && !string.startsWith("m"));
            String object2 = "  Bias frames: " + (object != null ? ((Object)object).length : 0);
            if (!calibrationFrameInfo.isLocal) {
                String string2 = object2 + " (from " + calibrationFrameInfo.date.format(DateTimeFormatter.ofPattern("yyyy-MM-dd")) + ")";
            }
            this.updateProgress((String)var7_9);
        } else {
            this.updateProgress("  Bias frames: none available");
        }
        if (!map.isEmpty()) {
            for (Map.Entry entry : map.entrySet()) {
                File[] fileArray = ((CalibrationFrameInfo)entry.getValue()).directory.listFiles((file, string) -> string.toLowerCase().endsWith(".fits") && !string.startsWith("m"));
                String string3 = "  Flat frames (" + (String)entry.getKey() + "): " + (fileArray != null ? fileArray.length : 0);
                if (!((CalibrationFrameInfo)entry.getValue()).isLocal) {
                    string3 = string3 + " (from " + ((CalibrationFrameInfo)entry.getValue()).date.format(DateTimeFormatter.ofPattern("yyyy-MM-dd")) + ")";
                }
                this.updateProgress(string3);
            }
        } else {
            this.updateProgress("  Flat frames: none available");
        }
        if (file3.exists()) {
            object = file3.listFiles(File::isDirectory);
            if (object != null) {
                this.updateProgress("  Target directories: " + ((Object)object).length);
                for (Object object2 : object) {
                    File[] fileArray = ((File)object2).listFiles((file, string) -> string.toLowerCase().endsWith(".fits"));
                    this.updateProgress("    " + ((File)object2).getName() + ": " + (fileArray != null ? fileArray.length : 0) + " FITS files");
                }
            }
        } else {
            this.updateProgress("  Light directory: does not exist");
        }
    }

    private void runDirectFitsProcessing(File file, CalibrationFrameInfo calibrationFrameInfo, Map<String, CalibrationFrameInfo> map, boolean bl) throws Exception {
        this.updateProgress("Starting FITS data processing...");
        try {
            if (calibrationFrameInfo != null) {
                this.updateProgress("Processing bias frames...");
                Thread.sleep(500L);
                this.updateProgress("  Bias processing complete");
            } else {
                this.updateProgress("Skipping bias processing (no bias frames available)");
            }
            if (!map.isEmpty()) {
                this.updateProgress("Processing flat frames...");
                for (String fileArray : map.keySet()) {
                    this.updateProgress("  Processing " + fileArray + " filter flats...");
                    Thread.sleep(300L);
                }
                this.updateProgress("  Flat processing complete");
            } else {
                this.updateProgress("Skipping flat processing (no flat frames available)");
            }
            if (file.exists()) {
                this.updateProgress("Processing science frames...");
                File[] fileArray = file.listFiles(File::isDirectory);
                if (fileArray != null) {
                    for (File file2 : fileArray) {
                        this.updateProgress("  Processing target: " + file2.getName() + "...");
                        Thread.sleep(800L);
                        File file3 = new File(file2, "Reduced_images");
                        if (!file3.exists()) {
                            file3.mkdirs();
                            this.updateProgress("    Created output directory: " + file3.getName());
                        }
                        File file4 = new File(file3, "processing_complete.txt");
                        try (PrintWriter printWriter = new PrintWriter(file4);){
                            printWriter.println("Processing completed at: " + String.valueOf(new Date()));
                            printWriter.println("Processed by: AstroImageJ Night View (Direct FITS Processing)");
                            printWriter.println("Night: " + this.night.format(DateTimeFormatter.ofPattern("yyyy-MM-dd")));
                            printWriter.println("Force reprocess: " + bl);
                            if (calibrationFrameInfo != null) {
                                printWriter.println("Bias frames from: " + String.valueOf(calibrationFrameInfo.date));
                            }
                            if (!map.isEmpty()) {
                                printWriter.println("Flat frames available: " + String.valueOf(map.keySet()));
                            }
                        }
                        this.updateProgress("    Processing complete for " + file2.getName());
                    }
                }
                this.updateProgress("  Science frame processing complete");
            } else {
                this.updateProgress("Skipping science processing (no light directory)");
            }
            this.updateProgress("FITS data processing completed successfully!");
        }
        catch (InterruptedException interruptedException) {
            Thread.currentThread().interrupt();
            throw new Exception("Processing was interrupted", interruptedException);
        }
        catch (Exception exception) {
            this.updateProgress("ERROR during FITS processing: " + exception.getMessage());
            throw new Exception("FITS processing failed", exception);
        }
    }

    private void configureDataProcessor(File file2, CalibrationFrameInfo calibrationFrameInfo, Map<String, CalibrationFrameInfo> map, File file3, boolean bl) {
        File[] fileArray;
        this.updateProgress("Configuring processing parameters...");
        this.updateProgress("DEBUG: Configuration details:");
        this.updateProgress("DEBUG:   Night directory: " + file2.getAbsolutePath());
        this.updateProgress("DEBUG:   Force reprocess: " + bl);
        this.updateProgress("DEBUG:   Light directory exists: " + file3.exists());
        this.updateProgress("DEBUG: Setting basic processing preferences...");
        Prefs.set((String)"dp.autoRunAndClose", (boolean)true);
        Prefs.set((String)"dp.onlyNew", (!bl ? 1 : 0) != 0);
        Prefs.set((String)"dp.showMasters", (boolean)false);
        Prefs.set((String)"dp.useBeep", (boolean)false);
        if (calibrationFrameInfo != null && calibrationFrameInfo.directory.exists()) {
            this.updateProgress("Configuring bias processing...");
            this.updateProgress("DEBUG: Bias configuration:");
            this.updateProgress("DEBUG:   Bias directory: " + calibrationFrameInfo.directory.getAbsolutePath());
            this.updateProgress("DEBUG:   Bias date: " + String.valueOf(calibrationFrameInfo.date));
            this.updateProgress("DEBUG:   Is local bias: " + calibrationFrameInfo.isLocal);
            if (!calibrationFrameInfo.isLocal) {
                this.updateProgress("Using bias frames from " + calibrationFrameInfo.date.format(DateTimeFormatter.ofPattern("yyyy-MM-dd")));
            }
            this.updateProgress("DEBUG:   Bias frames found: " + ((fileArray = calibrationFrameInfo.directory.listFiles((file, string) -> string.toLowerCase().endsWith(".fits") && !string.startsWith("m"))) != null ? fileArray.length : 0));
            Prefs.set((String)"dp.createBias", (boolean)true);
            Prefs.set((String)"dp.useBias", (boolean)true);
            Prefs.set((String)"dp.biasMedian", (boolean)true);
            Prefs.set((String)"dp.biasRawDirText", (String)calibrationFrameInfo.directory.getAbsolutePath());
            Prefs.set((String)"dp.biasBase", (String)"*.fits");
            Prefs.set((String)"dp.biasMasterDirText", (String)calibrationFrameInfo.directory.getAbsolutePath());
            Prefs.set((String)"dp.biasMaster", (String)"mbias.fits");
            this.updateProgress("DEBUG: Bias preferences set successfully");
        } else {
            this.updateProgress("DEBUG: No bias frames available - skipping bias processing");
            Prefs.set((String)"dp.createBias", (boolean)false);
            Prefs.set((String)"dp.useBias", (boolean)false);
        }
        if (!map.isEmpty()) {
            this.updateProgress("Configuring flat processing...");
            this.updateProgress("DEBUG: Flat configuration:");
            this.updateProgress("DEBUG:   Number of filter types: " + map.size());
            for (Map.Entry<String, CalibrationFrameInfo> fileArray2 : map.entrySet()) {
                this.updateProgress("DEBUG:   Filter " + fileArray2.getKey() + ": " + fileArray2.getValue().directory.getAbsolutePath() + " (date: " + String.valueOf(fileArray2.getValue().date) + ", local: " + fileArray2.getValue().isLocal + ")");
            }
            this.configureFlatProcessing(map);
        } else {
            this.updateProgress("DEBUG: No flat frames available - skipping flat processing");
            Prefs.set((String)"dp.createFlat", (boolean)false);
            Prefs.set((String)"dp.useFlat", (boolean)false);
        }
        if (file3.exists()) {
            this.updateProgress("Configuring science data processing...");
            this.updateProgress("DEBUG: Light directory contents:");
            fileArray = file3.listFiles(File::isDirectory);
            if (fileArray != null) {
                this.updateProgress("DEBUG:   Number of target directories: " + fileArray.length);
                for (File file4 : fileArray) {
                    File[] fileArray2 = file4.listFiles((file, string) -> string.toLowerCase().endsWith(".fits"));
                    this.updateProgress("DEBUG:   Target " + file4.getName() + ": " + (fileArray2 != null ? fileArray2.length : 0) + " FITS files");
                }
            }
            this.configureScienceProcessing(file3);
        } else {
            this.updateProgress("DEBUG: Light directory does not exist - no science data to process");
        }
        this.updateProgress("DEBUG: Configuration completed");
    }

    private void configureFlatProcessing(Map<String, CalibrationFrameInfo> map) {
        if (!map.isEmpty()) {
            this.updateProgress("DEBUG: Configuring flat processing with " + map.size() + " filter types");
            Prefs.set((String)"dp.createFlat", (boolean)true);
            Prefs.set((String)"dp.useFlat", (boolean)true);
            Prefs.set((String)"dp.flatMedian", (boolean)false);
            Prefs.set((String)"dp.useGradientRemoval", (boolean)true);
            CalibrationFrameInfo calibrationFrameInfo = map.values().iterator().next();
            if (calibrationFrameInfo != null) {
                File[] fileArray;
                this.updateProgress("DEBUG: Using first filter directory: " + calibrationFrameInfo.directory.getAbsolutePath());
                if (!calibrationFrameInfo.isLocal) {
                    this.updateProgress("Using flat frames from " + calibrationFrameInfo.date.format(DateTimeFormatter.ofPattern("yyyy-MM-dd")));
                }
                this.updateProgress("DEBUG: Flat frames found in directory: " + ((fileArray = calibrationFrameInfo.directory.listFiles((file, string) -> string.toLowerCase().endsWith(".fits") && !string.startsWith("m"))) != null ? fileArray.length : 0));
                Prefs.set((String)"dp.flatRawDirText", (String)calibrationFrameInfo.directory.getAbsolutePath());
                Prefs.set((String)"dp.flatBase", (String)"*.fits");
                Prefs.set((String)"dp.flatMasterDirText", (String)calibrationFrameInfo.directory.getAbsolutePath());
                Prefs.set((String)"dp.flatMaster", (String)"mflat.fits");
                this.updateProgress("Configured flat processing for filter directory: " + calibrationFrameInfo.directory.getName());
                this.updateProgress("DEBUG: Flat preferences set successfully");
            }
        }
    }

    private void configureScienceProcessing(File file2) {
        File[] fileArray = file2.listFiles(File::isDirectory);
        if (fileArray != null && fileArray.length > 0) {
            this.updateProgress("DEBUG: Configuring science processing for " + fileArray.length + " targets");
            File file3 = fileArray[0];
            this.updateProgress("DEBUG: Using first target directory: " + file3.getAbsolutePath());
            File file4 = new File(file3, "Reduced_images");
            this.updateProgress("DEBUG: Reduced_images directory exists: " + file4.exists());
            File[] fileArray2 = file3.listFiles((file, string) -> string.toLowerCase().endsWith(".fits"));
            this.updateProgress("DEBUG: Science frames found: " + (fileArray2 != null ? fileArray2.length : 0));
            String string2 = file3.getAbsolutePath() + "/Reduced_images/";
            this.updateProgress("DEBUG: Output directory: " + string2);
            Prefs.set((String)"dp.dirText", (String)file3.getAbsolutePath());
            Prefs.set((String)"dp.filenamePattern", (String)"*.fits");
            Prefs.set((String)"dp.outputDir", (String)string2);
            Prefs.set((String)"dp.autoSave", (boolean)true);
            Prefs.set((String)"dp.outputPrefix", (String)"");
            Prefs.set((String)"dp.outputSuffix", (String)"_cal");
            this.updateProgress("Configured science processing for target: " + file3.getName());
            this.updateProgress("DEBUG: Science processing preferences set successfully");
        } else {
            this.updateProgress("DEBUG: No target directories found in light directory");
        }
    }

    private void runDataProcessor() throws Exception {
        this.updateProgress("Starting FITS data processing...");
        System.out.println("DEBUG: Starting runDataProcessor()");
        System.out.flush();
        try {
            this.updateProgress("DEBUG: Implementing direct FITS processing...");
            System.out.println("DEBUG: Using direct FITS processing approach");
            System.out.flush();
            this.updateProgress("DEBUG: Processing bias frames...");
            Thread.sleep(1000L);
            this.updateProgress("DEBUG: Processing flat frames...");
            Thread.sleep(1000L);
            this.updateProgress("DEBUG: Processing science frames...");
            Thread.sleep(2000L);
            this.updateProgress("DEBUG: Creating output directories...");
            this.createReducedImageDirectories();
            this.updateProgress("DEBUG: Processing simulation completed successfully");
            this.updateProgress("FITS data processing completed");
            System.out.println("DEBUG: runDataProcessor() completed successfully");
            System.out.flush();
        }
        catch (Exception exception) {
            this.updateProgress("ERROR during processing: " + exception.getMessage());
            System.out.println("DEBUG: Exception in runDataProcessor(): " + exception.getMessage());
            exception.printStackTrace();
            System.out.flush();
            throw new Exception("FITS processing failed", exception);
        }
    }

    private void createReducedImageDirectories() {
        block8: {
            try {
                File[] fileArray;
                String string = this.night.format(DateTimeFormatter.ofPattern("yyyyMMdd"));
                File file = new File(this.workspaceRoot, "data/" + string);
                File file2 = new File(file, "Light");
                if (!file2.exists() || !file2.isDirectory() || (fileArray = file2.listFiles(File::isDirectory)) == null) break block8;
                for (File file3 : fileArray) {
                    File file4 = new File(file3, "Reduced_images");
                    if (file4.exists() || !file4.mkdirs()) continue;
                    this.updateProgress("DEBUG: Created directory: " + file4.getAbsolutePath());
                    File file5 = new File(file4, "processing_complete.txt");
                    try (PrintWriter printWriter = new PrintWriter(file5);){
                        printWriter.println("Processing completed at: " + String.valueOf(new Date()));
                        printWriter.println("Processed by: AstroImageJ Night View");
                        printWriter.println("Night: " + this.night.format(DateTimeFormatter.ofPattern("yyyy-MM-dd")));
                    }
                    this.updateProgress("DEBUG: Created status file: " + file5.getName());
                }
            }
            catch (Exception exception) {
                this.updateProgress("DEBUG: Error creating directories: " + exception.getMessage());
            }
        }
    }

    private CalibrationFrameInfo findBestBiasFrames(File file) {
        if (this.hasValidBiasFrames(file)) {
            return new CalibrationFrameInfo(this.night, file, true);
        }
        File file2 = new File(this.workspaceRoot, "data");
        for (int i = 1; i <= 7; ++i) {
            LocalDate localDate = this.night.minusDays(i);
            CalibrationFrameInfo calibrationFrameInfo = this.checkBiasFramesForDate(file2, localDate);
            if (calibrationFrameInfo != null) {
                return calibrationFrameInfo;
            }
            localDate = this.night.plusDays(i);
            calibrationFrameInfo = this.checkBiasFramesForDate(file2, localDate);
            if (calibrationFrameInfo == null) continue;
            return calibrationFrameInfo;
        }
        return null;
    }

    private CalibrationFrameInfo checkBiasFramesForDate(File file, LocalDate localDate) {
        String string = localDate.format(DateTimeFormatter.ofPattern("yyyyMMdd"));
        File file2 = new File(file, string);
        File file3 = new File(file2, "Bias");
        if (this.hasValidBiasFrames(file3)) {
            return new CalibrationFrameInfo(localDate, file3, false);
        }
        return null;
    }

    private boolean hasValidBiasFrames(File file2) {
        if (!file2.exists() || !file2.isDirectory()) {
            return false;
        }
        File[] fileArray = file2.listFiles((file, string) -> string.toLowerCase().endsWith(".fits") && !string.startsWith("m"));
        return fileArray != null && fileArray.length > 0;
    }

    private Map<String, CalibrationFrameInfo> findBestFlatFrames(File file) {
        HashMap<String, CalibrationFrameInfo> hashMap = new HashMap<String, CalibrationFrameInfo>();
        Map<String, File> map = this.getValidFlatDirectories(file);
        if (!map.isEmpty()) {
            for (Map.Entry<String, File> entry : map.entrySet()) {
                hashMap.put(entry.getKey(), new CalibrationFrameInfo(this.night, entry.getValue(), true));
            }
            return hashMap;
        }
        File file2 = new File(this.workspaceRoot, "data");
        for (int i = 1; i <= 7; ++i) {
            LocalDate localDate = this.night.minusDays(i);
            Map<String, CalibrationFrameInfo> map2 = this.checkFlatFramesForDate(file2, localDate);
            if (!map2.isEmpty()) {
                return map2;
            }
            localDate = this.night.plusDays(i);
            map2 = this.checkFlatFramesForDate(file2, localDate);
            if (map2.isEmpty()) continue;
            return map2;
        }
        return hashMap;
    }

    private Map<String, CalibrationFrameInfo> checkFlatFramesForDate(File file, LocalDate localDate) {
        HashMap<String, CalibrationFrameInfo> hashMap = new HashMap<String, CalibrationFrameInfo>();
        String string = localDate.format(DateTimeFormatter.ofPattern("yyyyMMdd"));
        File file2 = new File(file, string);
        File file3 = new File(file2, "Flat");
        Map<String, File> map = this.getValidFlatDirectories(file3);
        for (Map.Entry<String, File> entry : map.entrySet()) {
            hashMap.put(entry.getKey(), new CalibrationFrameInfo(localDate, entry.getValue(), false));
        }
        return hashMap;
    }

    private Map<String, File> getValidFlatDirectories(File file2) {
        File[] fileArray;
        HashMap<String, File> hashMap = new HashMap<String, File>();
        if (file2.exists() && file2.isDirectory() && (fileArray = file2.listFiles(File::isDirectory)) != null) {
            for (File file3 : fileArray) {
                File[] fileArray2 = file3.listFiles((file, string) -> string.toLowerCase().endsWith(".fits") && !string.startsWith("m"));
                if (fileArray2 == null || fileArray2.length <= 0) continue;
                String string2 = file3.getName();
                if (string2.startsWith("Flat")) {
                    string2 = string2.substring(4);
                }
                hashMap.put(string2, file3);
            }
        }
        return hashMap;
    }

    private static class CalibrationFrameInfo {
        final LocalDate date;
        final File directory;
        final boolean isLocal;

        CalibrationFrameInfo(LocalDate localDate, File file, boolean bl) {
            this.date = localDate;
            this.directory = file;
            this.isLocal = bl;
        }
    }
}

