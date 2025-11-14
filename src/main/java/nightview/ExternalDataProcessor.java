/*
 * Decompiled with CFR 0.152.
 */
package nightview;

import java.io.BufferedReader;
import java.io.File;
import java.io.InputStreamReader;
import java.time.LocalDate;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
import java.util.HashMap;
import java.util.Map;
import java.util.function.Consumer;

public class ExternalDataProcessor {
    private File workspaceRoot;
    private LocalDate night;
    private Consumer<String> progressCallback;
    private static final String ASTROIMAGEJ_APP = "/Applications/AstroImageJ.app";
    private static final String ASTROIMAGEJ_JAVA = "/Applications/AstroImageJ.app/Contents/PlugIns/jre.jre/Contents/Home/bin/java";
    private static final String ASTROIMAGEJ_JAR = "/Applications/AstroImageJ.app/Contents/Resources/ij.jar";
    private static final String ASTROIMAGEJ_PLUGINS = "/Applications/AstroImageJ.app/Contents/Resources/plugins";

    public ExternalDataProcessor(File file, LocalDate localDate) {
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
        System.out.println("DEBUG: ExternalDataProcessor.processNight() called");
        System.out.flush();
        String string = this.night.format(DateTimeFormatter.ofPattern("yyyyMMdd"));
        File file = new File(this.workspaceRoot, "data/" + string);
        System.out.println("DEBUG: Night directory: " + file.getAbsolutePath());
        System.out.flush();
        if (!file.exists()) {
            throw new RuntimeException("Night directory not found: " + file.getAbsolutePath());
        }
        this.updateProgress("Starting external AstroImageJ processing...");
        try {
            System.out.println("DEBUG: Checking AstroImageJ at: /Applications/AstroImageJ.app");
            System.out.flush();
            this.updateProgress("Checking AstroImageJ installation...");
            if (!new File(ASTROIMAGEJ_APP).exists()) {
                throw new Exception("AstroImageJ not found at: /Applications/AstroImageJ.app");
            }
            System.out.println("DEBUG: AstroImageJ found");
            System.out.flush();
            System.out.println("DEBUG: Starting data analysis...");
            System.out.flush();
            this.updateProgress("Analyzing data structure...");
            DataAnalysis dataAnalysis = this.analyzeNightData(file);
            System.out.println("DEBUG: Data analysis complete");
            System.out.flush();
            System.out.println("DEBUG: Generating script...");
            System.out.flush();
            this.updateProgress("Generating processing script...");
            String string2 = this.generateProcessingScript(dataAnalysis, bl);
            System.out.println("DEBUG: Script generated at: " + string2);
            System.out.flush();
            System.out.println("DEBUG: Executing AstroImageJ...");
            System.out.flush();
            this.updateProgress("Executing AstroImageJ Data Processor...");
            this.executeAstroImageJScript(string2);
            System.out.println("DEBUG: AstroImageJ execution complete");
            System.out.flush();
            this.updateProgress("External processing completed successfully!");
        }
        catch (Exception exception) {
            System.out.println("DEBUG: Exception caught: " + exception.getMessage());
            exception.printStackTrace();
            System.out.flush();
            String string3 = "External processing failed: " + exception.getMessage();
            this.updateProgress(string3);
            throw exception;
        }
    }

    private DataAnalysis analyzeNightData(File file) throws Exception {
        File[] fileArray;
        Object object;
        Object object2;
        DataAnalysis dataAnalysis = new DataAnalysis();
        dataAnalysis.nightDir = file;
        File file2 = new File(file, "Bias");
        if (this.hasValidBiasFrames(file2)) {
            dataAnalysis.biasDir = file2;
            dataAnalysis.biasDate = this.night;
            this.updateProgress("Local bias frames found");
        } else {
            object2 = this.findNearbyBiasFrames();
            if (object2 != null) {
                dataAnalysis.biasDir = ((BiasInfo)object2).directory;
                dataAnalysis.biasDate = ((BiasInfo)object2).date;
                this.updateProgress("Using bias frames from " + String.valueOf(((BiasInfo)object2).date));
            } else {
                this.updateProgress("No bias frames available");
            }
        }
        object2 = new File(file, "Flat");
        if (((File)object2).exists() && ((File)object2).isDirectory() && (object = ((File)object2).listFiles(File::isDirectory)) != null && ((Object)object).length > 0) {
            fileArray = object;
            int n = fileArray.length;
            for (int i = 0; i < n; ++i) {
                File file3 = fileArray[i];
                if (!this.hasValidFrames(file3)) continue;
                String object3 = this.extractFilterName(file3.getName());
                dataAnalysis.flatDirs.put(object3, file3);
            }
            if (!dataAnalysis.flatDirs.isEmpty()) {
                this.updateProgress("Local flat frames found: " + String.valueOf(dataAnalysis.flatDirs.keySet()));
            }
        }
        if (dataAnalysis.flatDirs.isEmpty()) {
            object = this.findNearbyFlatFrames();
            fileArray = object.entrySet().iterator();
            while (fileArray.hasNext()) {
                Map.Entry entry = (Map.Entry)fileArray.next();
                dataAnalysis.flatDirs.put((String)entry.getKey(), ((FlatInfo)entry.getValue()).directory);
            }
            if (!dataAnalysis.flatDirs.isEmpty()) {
                this.updateProgress("Using nearby flat frames: " + String.valueOf(dataAnalysis.flatDirs.keySet()));
            } else {
                this.updateProgress("No flat frames available");
            }
        }
        if (((File)(object = new File(file, "Light"))).exists() && ((File)object).isDirectory() && (fileArray = ((File)object).listFiles(File::isDirectory)) != null) {
            for (File file3 : fileArray) {
                if (!this.hasValidFrames(file3)) continue;
                dataAnalysis.targetDirs.put(file3.getName(), file3);
            }
        }
        if (dataAnalysis.targetDirs.isEmpty()) {
            throw new Exception("No science target directories found with FITS files");
        }
        this.updateProgress("Found targets: " + String.valueOf(dataAnalysis.targetDirs.keySet()));
        return dataAnalysis;
    }

    private String generateProcessingScript(DataAnalysis dataAnalysis, boolean bl) throws Exception {
        return "PYTHON_PROCESSING";
    }

    private void executeAstroImageJScript(String string) throws Exception {
        this.updateProgress("Starting Python FITS processing...");
        String string2 = LocalDateTime.now().format(DateTimeFormatter.ofPattern("yyyyMMdd_HHmmss"));
        File file = new File(this.workspaceRoot, "processing_log_" + string2 + ".txt");
        String string3 = this.night.format(DateTimeFormatter.ofPattern("yyyyMMdd"));
        File file2 = new File(this.workspaceRoot, "data/" + string3);
        DataAnalysis dataAnalysis = this.analyzeNightData(file2);
        int n = 0;
        for (Map.Entry<String, File> entry : dataAnalysis.targetDirs.entrySet()) {
            String string4 = entry.getKey();
            File file3 = entry.getValue();
            this.updateProgress("Processing target: " + string4);
            ProcessBuilder processBuilder = new ProcessBuilder("python3", new File(this.workspaceRoot, "process_fits.py").getAbsolutePath(), file3.getAbsolutePath(), dataAnalysis.biasDir != null ? dataAnalysis.biasDir.getAbsolutePath() : "", !dataAnalysis.flatDirs.isEmpty() ? dataAnalysis.flatDirs.values().iterator().next().getAbsolutePath() : "");
            processBuilder.directory(this.workspaceRoot);
            processBuilder.redirectErrorStream(true);
            this.updateProgress("Executing Python processor for " + string4);
            Process process = processBuilder.start();
            try (BufferedReader bufferedReader = new BufferedReader(new InputStreamReader(process.getInputStream()));){
                String string5;
                while ((string5 = bufferedReader.readLine()) != null) {
                    this.updateProgress(string5);
                }
            }
            int n2 = process.waitFor();
            if (n2 != 0) {
                throw new Exception("Python processing failed for " + string4 + " with exit code: " + n2);
            }
            ++n;
        }
        this.updateProgress("Processing completed - " + n + " targets processed");
        this.updateProgress("Log file: " + file.getName());
    }

    private boolean hasValidBiasFrames(File file2) {
        if (!file2.exists() || !file2.isDirectory()) {
            return false;
        }
        File[] fileArray = file2.listFiles((file, string) -> string.toLowerCase().endsWith(".fits") && !string.startsWith("m"));
        if (fileArray != null && fileArray.length > 0) {
            return true;
        }
        File[] fileArray2 = file2.listFiles((file, string) -> string.toLowerCase().endsWith(".fits") && string.startsWith("m"));
        return fileArray2 != null && fileArray2.length > 0;
    }

    private boolean hasValidFrames(File file2) {
        if (!file2.exists() || !file2.isDirectory()) {
            return false;
        }
        File[] fileArray = file2.listFiles((file, string) -> string.toLowerCase().endsWith(".fits") && !string.startsWith("m"));
        return fileArray != null && fileArray.length > 0;
    }

    private String extractFilterName(String string) {
        if (string.startsWith("Flat")) {
            return string.substring(4);
        }
        return string;
    }

    private BiasInfo findNearbyBiasFrames() {
        File file = new File(this.workspaceRoot, "data");
        for (int i = 1; i <= 7; ++i) {
            LocalDate localDate = this.night.minusDays(i);
            BiasInfo biasInfo = this.checkBiasFramesForDate(file, localDate);
            if (biasInfo != null) {
                return biasInfo;
            }
            localDate = this.night.plusDays(i);
            biasInfo = this.checkBiasFramesForDate(file, localDate);
            if (biasInfo == null) continue;
            return biasInfo;
        }
        return null;
    }

    private BiasInfo checkBiasFramesForDate(File file, LocalDate localDate) {
        String string = localDate.format(DateTimeFormatter.ofPattern("yyyyMMdd"));
        File file2 = new File(file, string);
        File file3 = new File(file2, "Bias");
        if (this.hasValidBiasFrames(file3)) {
            return new BiasInfo(localDate, file3);
        }
        return null;
    }

    private Map<String, FlatInfo> findNearbyFlatFrames() {
        File file = new File(this.workspaceRoot, "data");
        for (int i = 1; i <= 7; ++i) {
            LocalDate localDate = this.night.minusDays(i);
            Map<String, FlatInfo> map = this.checkFlatFramesForDate(file, localDate);
            if (!map.isEmpty()) {
                return map;
            }
            localDate = this.night.plusDays(i);
            map = this.checkFlatFramesForDate(file, localDate);
            if (map.isEmpty()) continue;
            return map;
        }
        return new HashMap<String, FlatInfo>();
    }

    private Map<String, FlatInfo> checkFlatFramesForDate(File file, LocalDate localDate) {
        File[] fileArray;
        HashMap<String, FlatInfo> hashMap = new HashMap<String, FlatInfo>();
        String string = localDate.format(DateTimeFormatter.ofPattern("yyyyMMdd"));
        File file2 = new File(file, string);
        File file3 = new File(file2, "Flat");
        if (file3.exists() && file3.isDirectory() && (fileArray = file3.listFiles(File::isDirectory)) != null) {
            for (File file4 : fileArray) {
                if (!this.hasValidFrames(file4)) continue;
                String string2 = this.extractFilterName(file4.getName());
                hashMap.put(string2, new FlatInfo(localDate, file4));
            }
        }
        return hashMap;
    }

    private static class DataAnalysis {
        File nightDir;
        File biasDir;
        LocalDate biasDate;
        Map<String, File> flatDirs = new HashMap<String, File>();
        Map<String, File> targetDirs = new HashMap<String, File>();

        private DataAnalysis() {
        }
    }

    private static class BiasInfo {
        LocalDate date;
        File directory;

        BiasInfo(LocalDate localDate, File file) {
            this.date = localDate;
            this.directory = file;
        }
    }

    private static class FlatInfo {
        LocalDate date;
        File directory;

        FlatInfo(LocalDate localDate, File file) {
            this.date = localDate;
            this.directory = file;
        }
    }
}

