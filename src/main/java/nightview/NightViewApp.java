/*
 * Decompiled with CFR 0.152.
 */
package nightview;

import java.awt.BorderLayout;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.Font;
import java.awt.Frame;
import java.io.File;
import java.time.LocalDate;
import java.time.format.DateTimeFormatter;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTabbedPane;
import javax.swing.JTextArea;
import javax.swing.SwingUtilities;
import javax.swing.SwingWorker;
import javax.swing.Timer;
import javax.swing.UIManager;
import nightview.ExternalDataProcessor;
import nightview.NightCalendarPanel;
import nightview.ProgressDialog;
import nightview.TargetPhotometryPanel;

public class NightViewApp
extends JFrame {
    private File workspaceRoot = new File(System.getProperty("user.dir"));
    private LocalDate selectedNight;
    private NightCalendarPanel calendarPanel;
    private JPanel nightDetailsPanel;
    private JLabel nightLabel;
    private JTextArea detailsArea;
    private JButton processButton;
    private JButton reprocessButton;

    public NightViewApp() {
        System.out.println("Workspace root: " + this.workspaceRoot.getAbsolutePath());
        this.initializeGUI();
        this.setSelectedNight(LocalDate.now());
    }

    private void initializeGUI() {
        this.setTitle("AstroImageJ Night View");
        this.setDefaultCloseOperation(3);
        this.setLayout(new BorderLayout());
        JTabbedPane jTabbedPane = new JTabbedPane();
        JPanel jPanel = new JPanel(new BorderLayout());
        this.calendarPanel = new NightCalendarPanel(this);
        this.createNightDetailsPanel();
        jPanel.add((Component)this.calendarPanel, "West");
        jPanel.add((Component)this.nightDetailsPanel, "Center");
        jTabbedPane.addTab("Night View", jPanel);
        TargetPhotometryPanel targetPhotometryPanel = new TargetPhotometryPanel(this.workspaceRoot);
        jTabbedPane.addTab("Target Photometry", targetPhotometryPanel);
        this.add((Component)jTabbedPane, "Center");
        this.pack();
        this.setLocationRelativeTo(null);
    }

    private void createNightDetailsPanel() {
        this.nightDetailsPanel = new JPanel(new BorderLayout());
        this.nightDetailsPanel.setBorder(BorderFactory.createTitledBorder("Night Details"));
        this.nightDetailsPanel.setPreferredSize(new Dimension(400, 500));
        this.nightLabel = new JLabel("No night selected", 0);
        this.nightLabel.setFont(this.nightLabel.getFont().deriveFont(1, 16.0f));
        this.nightLabel.setBorder(BorderFactory.createEmptyBorder(10, 10, 10, 10));
        this.detailsArea = new JTextArea(20, 30);
        this.detailsArea.setEditable(false);
        this.detailsArea.setFont(new Font("Monospaced", 0, 12));
        JScrollPane jScrollPane = new JScrollPane(this.detailsArea);
        JPanel jPanel = new JPanel(new FlowLayout());
        this.processButton = new JButton("Process Night");
        this.reprocessButton = new JButton("Reprocess Night");
        this.processButton.addActionListener(actionEvent -> this.processNight());
        this.reprocessButton.addActionListener(actionEvent -> this.reprocessNight());
        jPanel.add(this.processButton);
        jPanel.add(this.reprocessButton);
        this.nightDetailsPanel.add((Component)this.nightLabel, "North");
        this.nightDetailsPanel.add((Component)jScrollPane, "Center");
        this.nightDetailsPanel.add((Component)jPanel, "South");
    }

    public void setSelectedNight(LocalDate localDate) {
        this.selectedNight = localDate;
        this.updateNightDetails();
    }

    public void onNightSelected(LocalDate localDate) {
        this.setSelectedNight(localDate);
    }

    private void updateNightDetails() {
        if (this.selectedNight == null) {
            this.nightLabel.setText("No night selected");
            this.detailsArea.setText("");
            this.processButton.setEnabled(false);
            this.reprocessButton.setEnabled(false);
            return;
        }
        this.nightLabel.setText("Night: " + this.selectedNight.format(DateTimeFormatter.ofPattern("yyyy-MM-dd")));
        this.calendarPanel.setAvailableNights(this.getAllAvailableNights());
        this.calendarPanel.setProcessedNights(this.getAllProcessedNights());
        String string = this.analyzeNightData(this.selectedNight);
        this.detailsArea.setText(string);
        String string2 = this.selectedNight.format(DateTimeFormatter.ofPattern("yyyyMMdd"));
        File file = new File(this.workspaceRoot, "data/" + string2);
        boolean bl = this.isNightFullyProcessed(file);
        this.processButton.setEnabled(!bl);
        this.reprocessButton.setEnabled(true);
    }

    private Set<LocalDate> getAllAvailableNights() {
        HashSet<LocalDate> hashSet = new HashSet<LocalDate>();
        File file = new File(this.workspaceRoot, "data");
        if (file.exists() && file.isDirectory()) {
            for (File file2 : file.listFiles()) {
                if (!file2.isDirectory()) continue;
                try {
                    File[] fileArray;
                    File file3;
                    File file4;
                    LocalDate localDate = LocalDate.parse(file2.getName(), DateTimeFormatter.ofPattern("yyyyMMdd"));
                    boolean bl = false;
                    File file5 = new File(file2, "Bias");
                    if (file5.exists() && this.hasAnyFitsFiles(file5)) {
                        bl = true;
                    }
                    if ((file4 = new File(file2, "Flat")).exists() && this.hasAnyFitsFiles(file4)) {
                        bl = true;
                    }
                    if ((file3 = new File(file2, "Light")).exists() && (fileArray = file3.listFiles(File::isDirectory)) != null) {
                        for (File file6 : fileArray) {
                            if (!this.hasAnyFitsFiles(file6)) continue;
                            bl = true;
                            break;
                        }
                    }
                    if (!bl) continue;
                    hashSet.add(localDate);
                }
                catch (Exception exception) {
                    // empty catch block
                }
            }
        }
        return hashSet;
    }

    private Set<LocalDate> getAllProcessedNights() {
        HashSet<LocalDate> hashSet = new HashSet<LocalDate>();
        File file = new File(this.workspaceRoot, "data");
        if (file.exists() && file.isDirectory()) {
            for (File file2 : file.listFiles()) {
                if (!file2.isDirectory()) continue;
                try {
                    LocalDate localDate = LocalDate.parse(file2.getName(), DateTimeFormatter.ofPattern("yyyyMMdd"));
                    if (!this.isNightFullyProcessed(file2)) continue;
                    hashSet.add(localDate);
                }
                catch (Exception exception) {
                    // empty catch block
                }
            }
        }
        return hashSet;
    }

    private boolean hasAnyFitsFiles(File file2) {
        if (!file2.exists() || !file2.isDirectory()) {
            return false;
        }
        File[] fileArray = file2.listFiles((file, string) -> string.toLowerCase().endsWith(".fits"));
        if (fileArray != null && fileArray.length > 0) {
            return true;
        }
        File[] fileArray2 = file2.listFiles(File::isDirectory);
        if (fileArray2 != null) {
            for (File file3 : fileArray2) {
                File[] fileArray3 = file3.listFiles((file, string) -> string.toLowerCase().endsWith(".fits"));
                if (fileArray3 == null || fileArray3.length <= 0) continue;
                return true;
            }
        }
        return false;
    }

    private Set<LocalDate> getProcessedNights() {
        HashSet<LocalDate> hashSet = new HashSet<LocalDate>();
        File file = new File(this.workspaceRoot, "data");
        if (file.exists() && file.isDirectory()) {
            block2: for (File file2 : file.listFiles()) {
                if (!file2.isDirectory()) continue;
                try {
                    File[] fileArray;
                    LocalDate localDate = LocalDate.parse(file2.getName(), DateTimeFormatter.ofPattern("yyyyMMdd"));
                    File file3 = new File(file2, "Light");
                    if (!file3.exists() || (fileArray = file3.listFiles(File::isDirectory)) == null) continue;
                    for (File file4 : fileArray) {
                        File file5 = new File(file4, "Reduced_images");
                        if (!file5.exists()) continue;
                        hashSet.add(localDate);
                        continue block2;
                    }
                }
                catch (Exception exception) {
                    // empty catch block
                }
            }
        }
        return hashSet;
    }

    private String analyzeNightData(LocalDate localDate) {
        StringBuilder stringBuilder = new StringBuilder();
        String string = localDate.format(DateTimeFormatter.ofPattern("yyyyMMdd"));
        File file = new File(this.workspaceRoot, "data/" + string);
        stringBuilder.append("Data Analysis for ").append(string).append("\n");
        stringBuilder.append("=====================================\n\n");
        if (!file.exists()) {
            stringBuilder.append("Night directory not found: ").append(file.getAbsolutePath()).append("\n");
            return stringBuilder.toString();
        }
        stringBuilder.append("Night directory: ").append(file.getAbsolutePath()).append("\n\n");
        this.analyzeBiasFrames(file, stringBuilder, localDate);
        this.analyzeFlatFrames(file, stringBuilder, localDate);
        this.analyzeLightFrames(file, stringBuilder);
        return stringBuilder.toString();
    }

    private void analyzeBiasFrames(File file2, StringBuilder stringBuilder, LocalDate localDate) {
        File file3 = new File(file2, "Bias");
        stringBuilder.append("BIAS FRAMES:\n");
        if (!file3.exists()) {
            BiasInfo biasInfo = this.findNearbyBiasFrames(localDate);
            if (biasInfo != null) {
                stringBuilder.append("  No local bias frames found\n");
                stringBuilder.append("  Using bias frames from: ").append(biasInfo.date.format(DateTimeFormatter.ofPattern("yyyy-MM-dd"))).append("\n");
                stringBuilder.append("  Location: ").append(biasInfo.directory.getAbsolutePath()).append("\n");
                File[] fileArray = biasInfo.directory.listFiles((file, string) -> string.toLowerCase().endsWith(".fits") && !string.startsWith("m"));
                stringBuilder.append("  Bias frames: ").append(fileArray != null ? fileArray.length : 0).append("\n");
                File[] fileArray2 = biasInfo.directory.listFiles((file, string) -> string.toLowerCase().equals("mbias.fits"));
                stringBuilder.append("  Master bias: ").append(fileArray2 != null ? fileArray2.length : 0).append("\n");
            } else {
                stringBuilder.append("  No bias frames found (searched +/- 7 days)\n");
            }
        } else {
            File[] fileArray = file3.listFiles((file, string) -> string.toLowerCase().endsWith(".fits") && !string.startsWith("m"));
            stringBuilder.append("  Local bias frames available\n");
            stringBuilder.append("  Location: ").append(file3.getAbsolutePath()).append("\n");
            stringBuilder.append("  Bias frames: ").append(fileArray != null ? fileArray.length : 0).append("\n");
            File[] fileArray3 = file3.listFiles((file, string) -> string.toLowerCase().equals("mbias.fits"));
            stringBuilder.append("  Master bias: ").append(fileArray3 != null ? fileArray3.length : 0).append("\n");
        }
        stringBuilder.append("\n");
    }

    private void analyzeFlatFrames(File file2, StringBuilder stringBuilder, LocalDate localDate) {
        File file3 = new File(file2, "Flat");
        stringBuilder.append("FLAT FRAMES:\n");
        if (!file3.exists()) {
            Map<String, FlatInfo> map = this.findNearbyFlatFrames(localDate);
            if (!map.isEmpty()) {
                stringBuilder.append("  No local flat frames found\n");
                for (Map.Entry<String, FlatInfo> entry : map.entrySet()) {
                    FlatInfo flatInfo = entry.getValue();
                    stringBuilder.append("  Using ").append(entry.getKey()).append(" flats from: ").append(flatInfo.date.format(DateTimeFormatter.ofPattern("yyyy-MM-dd"))).append("\n");
                    stringBuilder.append("  Location: ").append(flatInfo.directory.getAbsolutePath()).append("\n");
                    File[] fileArray = flatInfo.directory.listFiles((file, string) -> string.toLowerCase().endsWith(".fits") && !string.startsWith("m"));
                    stringBuilder.append("  Flat frames: ").append(fileArray != null ? fileArray.length : 0).append("\n");
                    File[] fileArray2 = flatInfo.directory.listFiles((file, string) -> string.toLowerCase().equals("mflat.fits"));
                    stringBuilder.append("  Master flat: ").append(fileArray2 != null ? fileArray2.length : 0).append("\n");
                }
            } else {
                stringBuilder.append("  No flat frames found (searched +/- 7 days)\n");
            }
        } else {
            File[] fileArray = file3.listFiles(File::isDirectory);
            if (fileArray != null && fileArray.length > 0) {
                stringBuilder.append("  Local flat frames available\n");
                for (File file4 : fileArray) {
                    String string2 = file4.getName();
                    File[] fileArray3 = file4.listFiles((file, string) -> string.toLowerCase().endsWith(".fits") && !string.startsWith("m"));
                    stringBuilder.append("  ").append(string2).append(" frames: ").append(fileArray3 != null ? fileArray3.length : 0).append("\n");
                    File[] fileArray4 = file4.listFiles((file, string) -> string.toLowerCase().equals("mflat.fits"));
                    stringBuilder.append("  ").append(string2).append(" master: ").append(fileArray4 != null ? fileArray4.length : 0).append("\n");
                }
            } else {
                stringBuilder.append("  No filter directories found\n");
            }
        }
        stringBuilder.append("\n");
    }

    private void analyzeLightFrames(File file2, StringBuilder stringBuilder) {
        File file3 = new File(file2, "Light");
        stringBuilder.append("LIGHT FRAMES:\n");
        if (!file3.exists()) {
            stringBuilder.append("  Light directory not found\n");
            return;
        }
        File[] fileArray = file3.listFiles(File::isDirectory);
        if (fileArray == null || fileArray.length == 0) {
            stringBuilder.append("  No target directories found\n");
            return;
        }
        stringBuilder.append("  Target directories found: ").append(fileArray.length).append("\n");
        for (File file4 : fileArray) {
            stringBuilder.append("  ").append(file4.getName()).append(":\n");
            File[] fileArray2 = file4.listFiles((file, string) -> string.toLowerCase().endsWith(".fits"));
            int n = fileArray2 != null ? fileArray2.length : 0;
            stringBuilder.append("    FITS files: ").append(n).append("\n");
            File file5 = new File(file4, "Reduced_images");
            if (file5.exists()) {
                File[] fileArray3 = file5.listFiles((file, string) -> string.toLowerCase().endsWith(".fits"));
                int n2 = fileArray3 != null ? fileArray3.length : 0;
                stringBuilder.append("    Processed files: ").append(n2).append("\n");
                if (n > 0 && n2 >= n) {
                    stringBuilder.append("    Status: FULLY PROCESSED\n");
                    continue;
                }
                stringBuilder.append("    Status: ").append(n2).append("/").append(n).append(" processed\n");
                continue;
            }
            stringBuilder.append("    Not yet processed\n");
        }
        stringBuilder.append("\n");
    }

    private boolean isNightFullyProcessed(File file2) {
        File file3 = new File(file2, "Light");
        if (!file3.exists()) {
            return false;
        }
        File[] fileArray = file3.listFiles(File::isDirectory);
        if (fileArray == null || fileArray.length == 0) {
            return false;
        }
        for (File file4 : fileArray) {
            int n;
            int n2;
            File[] fileArray2 = file4.listFiles((file, string) -> string.toLowerCase().endsWith(".fits"));
            int n3 = n2 = fileArray2 != null ? fileArray2.length : 0;
            if (n2 == 0) continue;
            File file5 = new File(file4, "Reduced_images");
            if (!file5.exists()) {
                return false;
            }
            File[] fileArray3 = file5.listFiles((file, string) -> string.toLowerCase().endsWith(".fits"));
            int n4 = n = fileArray3 != null ? fileArray3.length : 0;
            if (n >= n2) continue;
            return false;
        }
        return true;
    }

    private BiasInfo findNearbyBiasFrames(LocalDate localDate) {
        File file = new File(this.workspaceRoot, "data");
        for (int i = 1; i <= 7; ++i) {
            LocalDate localDate2 = localDate.minusDays(i);
            BiasInfo biasInfo = this.checkBiasFramesForDate(file, localDate2);
            if (biasInfo != null) {
                return biasInfo;
            }
            localDate2 = localDate.plusDays(i);
            biasInfo = this.checkBiasFramesForDate(file, localDate2);
            if (biasInfo == null) continue;
            return biasInfo;
        }
        return null;
    }

    private BiasInfo checkBiasFramesForDate(File file2, LocalDate localDate) {
        String string2 = localDate.format(DateTimeFormatter.ofPattern("yyyyMMdd"));
        File file3 = new File(file2, string2);
        File file4 = new File(file3, "Bias");
        if (file4.exists() && file4.isDirectory()) {
            File[] fileArray = file4.listFiles((file, string) -> string.toLowerCase().endsWith(".fits") && !string.startsWith("m"));
            if (fileArray != null && fileArray.length > 0) {
                return new BiasInfo(localDate, file4);
            }
            File[] fileArray2 = file4.listFiles((file, string) -> string.toLowerCase().endsWith(".fits") && string.startsWith("m"));
            if (fileArray2 != null && fileArray2.length > 0) {
                return new BiasInfo(localDate, file4);
            }
        }
        return null;
    }

    private Map<String, FlatInfo> findNearbyFlatFrames(LocalDate localDate) {
        File file = new File(this.workspaceRoot, "data");
        for (int i = 1; i <= 7; ++i) {
            LocalDate localDate2 = localDate.minusDays(i);
            Map<String, FlatInfo> map = this.checkFlatFramesForDate(file, localDate2);
            if (!map.isEmpty()) {
                return map;
            }
            localDate2 = localDate.plusDays(i);
            map = this.checkFlatFramesForDate(file, localDate2);
            if (map.isEmpty()) continue;
            return map;
        }
        return new HashMap<String, FlatInfo>();
    }

    private Map<String, FlatInfo> checkFlatFramesForDate(File file2, LocalDate localDate) {
        File[] fileArray;
        HashMap<String, FlatInfo> hashMap = new HashMap<String, FlatInfo>();
        String string2 = localDate.format(DateTimeFormatter.ofPattern("yyyyMMdd"));
        File file3 = new File(file2, string2);
        File file4 = new File(file3, "Flat");
        if (file4.exists() && file4.isDirectory() && (fileArray = file4.listFiles(File::isDirectory)) != null) {
            for (File file5 : fileArray) {
                File[] fileArray2 = file5.listFiles((file, string) -> string.toLowerCase().endsWith(".fits") && !string.startsWith("m"));
                if (fileArray2 == null || fileArray2.length <= 0) continue;
                String string3 = file5.getName();
                if (string3.startsWith("Flat")) {
                    string3 = string3.substring(4);
                }
                hashMap.put(string3, new FlatInfo(localDate, file5));
            }
        }
        return hashMap;
    }

    private void processNight() {
        System.out.println("DEBUG: processNight() called");
        System.out.flush();
        if (this.selectedNight == null) {
            System.out.println("DEBUG: No night selected!");
            JOptionPane.showMessageDialog(this, "Please select a night first", "Error", 0);
            return;
        }
        System.out.println("DEBUG: Selected night: " + String.valueOf(this.selectedNight));
        System.out.flush();
        SwingUtilities.invokeLater(() -> {
            System.out.println("DEBUG: Inside invokeLater");
            System.out.flush();
            System.out.println("DEBUG: Creating ProgressDialog...");
            System.out.flush();
            final ProgressDialog progressDialog = new ProgressDialog((Frame)this, "Processing Night Data");
            System.out.println("DEBUG: ProgressDialog created, making visible...");
            System.out.flush();
            progressDialog.setVisible(true);
            System.out.println("DEBUG: ProgressDialog is visible");
            System.out.flush();
            Timer timer = new Timer(100, actionEvent -> {
                System.out.println("DEBUG: Creating SwingWorker...");
                System.out.flush();
                SwingWorker<Void, String> swingWorker = new SwingWorker<Void, String>(this){
                    final /* synthetic */ NightViewApp this$0;
                    {
                        this.this$0 = nightViewApp;
                    }

                    @Override
                    protected Void doInBackground() throws Exception {
                        System.out.println("DEBUG: SwingWorker.doInBackground() started");
                        System.out.flush();
                        try {
                            System.out.println("DEBUG: Creating ExternalDataProcessor...");
                            System.out.flush();
                            ExternalDataProcessor externalDataProcessor = new ExternalDataProcessor(this.this$0.workspaceRoot, this.this$0.selectedNight);
                            System.out.println("DEBUG: Setting progress callback...");
                            System.out.flush();
                            1 var2_3 = this;
                            externalDataProcessor.setProgressCallback(string -> var2_3.publish(string));
                            System.out.println("DEBUG: Calling processor.processNight()...");
                            System.out.flush();
                            externalDataProcessor.processNight(false);
                            System.out.println("DEBUG: processor.processNight() completed");
                            System.out.flush();
                        }
                        catch (Exception exception) {
                            System.out.println("DEBUG: Exception in doInBackground: " + exception.getMessage());
                            exception.printStackTrace();
                            System.out.flush();
                            throw exception;
                        }
                        return null;
                    }

                    @Override
                    protected void process(List<String> list) {
                        System.out.println("DEBUG: process() called with " + list.size() + " messages");
                        System.out.flush();
                        for (String string : list) {
                            progressDialog.updateProgress(string);
                        }
                    }

                    @Override
                    protected void done() {
                        System.out.println("DEBUG: SwingWorker.done() called");
                        System.out.flush();
                        progressDialog.setVisible(false);
                        try {
                            this.get();
                            System.out.println("DEBUG: Processing successful");
                            System.out.flush();
                            JOptionPane.showMessageDialog(this.this$0, "Processing completed successfully!", "Success", 1);
                            this.this$0.updateNightDetails();
                        }
                        catch (Exception exception) {
                            System.out.println("DEBUG: Processing failed with exception");
                            exception.printStackTrace();
                            System.out.flush();
                            JOptionPane.showMessageDialog(this.this$0, "Processing failed: " + exception.getMessage(), "Error", 0);
                        }
                    }
                };
                System.out.println("DEBUG: Executing SwingWorker...");
                System.out.flush();
                swingWorker.execute();
                System.out.println("DEBUG: SwingWorker.execute() called");
                System.out.flush();
            });
            timer.setRepeats(false);
            timer.start();
        });
    }

    private void reprocessNight() {
        SwingUtilities.invokeLater(() -> {
            final ProgressDialog progressDialog = new ProgressDialog((Frame)this, "Reprocessing Night Data");
            progressDialog.setVisible(true);
            SwingWorker<Void, String> swingWorker = new SwingWorker<Void, String>(this){
                final /* synthetic */ NightViewApp this$0;
                {
                    this.this$0 = nightViewApp;
                }

                @Override
                protected Void doInBackground() throws Exception {
                    ExternalDataProcessor externalDataProcessor = new ExternalDataProcessor(this.this$0.workspaceRoot, this.this$0.selectedNight);
                    2 var2_2 = this;
                    externalDataProcessor.setProgressCallback(string -> var2_2.publish(string));
                    externalDataProcessor.processNight(true);
                    return null;
                }

                @Override
                protected void process(List<String> list) {
                    for (String string : list) {
                        progressDialog.updateProgress(string);
                    }
                }

                @Override
                protected void done() {
                    progressDialog.setVisible(false);
                    try {
                        this.get();
                        JOptionPane.showMessageDialog(this.this$0, "Reprocessing completed successfully!", "Success", 1);
                        this.this$0.updateNightDetails();
                    }
                    catch (Exception exception) {
                        exception.printStackTrace();
                        JOptionPane.showMessageDialog(this.this$0, "Reprocessing failed: " + exception.getMessage(), "Error", 0);
                    }
                }
            };
            swingWorker.execute();
        });
    }

    public static void main(String[] stringArray) {
        SwingUtilities.invokeLater(() -> {
            try {
                for (UIManager.LookAndFeelInfo lookAndFeelInfo : UIManager.getInstalledLookAndFeels()) {
                    if (!"Nimbus".equals(lookAndFeelInfo.getName())) continue;
                    UIManager.setLookAndFeel(lookAndFeelInfo.getClassName());
                    break;
                }
            }
            catch (Exception exception) {
                exception.printStackTrace();
            }
            new NightViewApp().setVisible(true);
        });
    }

    private static class BiasInfo {
        final LocalDate date;
        final File directory;

        BiasInfo(LocalDate localDate, File file) {
            this.date = localDate;
            this.directory = file;
        }
    }

    private static class FlatInfo {
        final LocalDate date;
        final File directory;

        FlatInfo(LocalDate localDate, File file) {
            this.date = localDate;
            this.directory = file;
        }
    }
}

