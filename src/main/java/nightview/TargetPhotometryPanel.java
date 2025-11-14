/*
 * Decompiled with CFR 0.152.
 * 
 * Could not load the following classes:
 *  nom.tam.fits.Fits
 *  nom.tam.fits.Header
 *  nom.tam.fits.ImageData
 *  nom.tam.fits.ImageHDU
 */
package nightview;

import java.awt.BasicStroke;
import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Frame;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.GridLayout;
import java.awt.RenderingHints;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.MouseMotionAdapter;
import java.awt.image.BufferedImage;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.net.HttpURLConnection;
import java.net.URL;
import java.net.URLConnection;
import java.net.URLEncoder;
import java.time.LocalDate;
import java.time.format.DateTimeFormatter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSpinner;
import javax.swing.JSplitPane;
import javax.swing.JTabbedPane;
import javax.swing.JTable;
import javax.swing.JTextArea;
import javax.swing.SpinnerNumberModel;
import javax.swing.SwingUtilities;
import javax.swing.SwingWorker;
import javax.swing.table.DefaultTableCellRenderer;
import javax.swing.table.DefaultTableModel;
import nom.tam.fits.Fits;
import nom.tam.fits.Header;
import nom.tam.fits.ImageData;
import nom.tam.fits.ImageHDU;

public class TargetPhotometryPanel
extends JPanel {
    private File workspaceRoot;
    private JTable targetTable;
    private DefaultTableModel tableModel;
    private ImageDisplayPanel imagePanel;
    private JSpinner fovSpinner;
    private JSpinner magLimitSpinner;
    private JLabel statusLabel;
    private JPanel radialProfilePanel;
    private File currentFitsFile;
    private File currentTargetDir;
    private Header currentFitsHeader;
    private BufferedImage currentImage;
    private int[][] currentImageData;
    private List<ReferenceStar> referenceStars;
    private TargetInfo targetInfo;
    private List<TargetInfo> variableStarsInField;
    private double currentMagLimit = 14.5;
    private double currentFOV = 60.0;
    private double minVarDeltaV = 0.06;
    private StarLocation radialProfileStar = null;
    private RadialProfile currentRadialProfile = null;
    private List<PhotometryResult> photometryResults = null;
    private ReferenceStar checkStar = null;
    private JPanel photometryResultsPanel = null;
    private JPanel lightcurvePanel = null;
    private JComboBox<String> variableSelector = null;
    private Object selectedStar = null;
    private double selectedStarRA = 0.0;
    private double selectedStarDec = 0.0;
    private double zoomLevel = 1.0;
    private int panOffsetX = 0;
    private int panOffsetY = 0;
    private int dragStartX = 0;
    private int dragStartY = 0;
    private boolean isDragging = false;
    private JTextArea detailsArea;
    private JButton resetZoomButton;
    private JSpinner varDeltaVSpinner;
    private boolean showApertures = false;
    private JLabel mouseInfoLabel;
    private static final String STATUS_NOT_DONE = "Not Done";
    private static final String STATUS_DONE = "Done";
    private static final String STATUS_BAD_OBS = "Bad Obs";
    private static final String STATUS_FILE = "photometry_status.txt";

    public TargetPhotometryPanel(File file) {
        this.workspaceRoot = file;
        this.setLayout(new BorderLayout(10, 10));
        this.setBorder(BorderFactory.createEmptyBorder(10, 10, 10, 10));
        JSplitPane jSplitPane = new JSplitPane(1);
        jSplitPane.setDividerLocation(400);
        JPanel jPanel = this.createTargetListPanel();
        jSplitPane.setLeftComponent(jPanel);
        JPanel jPanel2 = this.createImageViewerPanel();
        jSplitPane.setRightComponent(jPanel2);
        this.add((Component)jSplitPane, "Center");
        this.refreshTargetList();
    }

    private JPanel createTargetListPanel() {
        JPanel jPanel = new JPanel(new BorderLayout(5, 5));
        jPanel.setBorder(BorderFactory.createTitledBorder("Target-Night List"));
        Object[] objectArray = new String[]{"Target", "Night", "Status", "Images"};
        this.tableModel = new DefaultTableModel(objectArray, 0){

            @Override
            public boolean isCellEditable(int n, int n2) {
                return false;
            }
        };
        this.targetTable = new JTable(this.tableModel);
        this.targetTable.setSelectionMode(0);
        this.targetTable.getSelectionModel().addListSelectionListener(listSelectionEvent -> {
            if (!listSelectionEvent.getValueIsAdjusting()) {
                this.onTargetSelected();
            }
        });
        this.targetTable.getColumnModel().getColumn(0).setPreferredWidth(120);
        this.targetTable.getColumnModel().getColumn(1).setPreferredWidth(100);
        this.targetTable.getColumnModel().getColumn(2).setPreferredWidth(80);
        this.targetTable.getColumnModel().getColumn(3).setPreferredWidth(60);
        this.targetTable.getColumnModel().getColumn(2).setCellRenderer(new StatusCellRenderer());
        JScrollPane jScrollPane = new JScrollPane(this.targetTable);
        jPanel.add((Component)jScrollPane, "Center");
        JButton jButton = new JButton("Refresh");
        jButton.addActionListener(actionEvent -> this.refreshTargetList());
        jPanel.add((Component)jButton, "South");
        return jPanel;
    }

    private JPanel createImageViewerPanel() {
        JPanel jPanel = new JPanel(new BorderLayout(5, 5));
        jPanel.setBorder(BorderFactory.createTitledBorder("Target Photometry"));
        JTabbedPane jTabbedPane = new JTabbedPane();
        JPanel jPanel2 = this.createImageViewerTab();
        jTabbedPane.addTab("Image Viewer", jPanel2);
        this.radialProfilePanel = this.createRadialProfileTab();
        jTabbedPane.addTab("Aperture Photometry Settings", this.radialProfilePanel);
        JPanel jPanel3 = this.createPhotometryResultsTab();
        jTabbedPane.addTab("Photometry Results", jPanel3);
        JPanel jPanel4 = this.createLightcurvesTab();
        jTabbedPane.addTab("Lightcurves", jPanel4);
        jPanel.add((Component)jTabbedPane, "Center");
        return jPanel;
    }

    private JPanel createImageViewerTab() {
        JPanel jPanel = new JPanel(new BorderLayout(5, 5));
        this.imagePanel = new ImageDisplayPanel();
        this.imagePanel.setBackground(Color.BLACK);
        this.imagePanel.setPreferredSize(new Dimension(600, 600));
        jPanel.add((Component)this.imagePanel, "Center");
        JPanel jPanel2 = new JPanel(new BorderLayout());
        JPanel jPanel3 = new JPanel(new BorderLayout(5, 5));
        JPanel jPanel4 = new JPanel(new FlowLayout(0));
        this.resetZoomButton = new JButton("Reset Zoom");
        this.resetZoomButton.addActionListener(actionEvent -> this.resetZoom());
        jPanel4.add(this.resetZoomButton);
        JButton jButton = new JButton("Do Phot All");
        jButton.addActionListener(actionEvent -> this.performPhotometryAll());
        jPanel4.add(jButton);
        JButton jButton2 = new JButton("Clear Results");
        jButton2.addActionListener(actionEvent -> {
            this.photometryResults = null;
            this.updatePhotometryResultsTable();
            this.statusLabel.setText("Results cleared");
            this.statusLabel.setForeground(Color.BLUE);
        });
        jPanel4.add(jButton2);
        this.statusLabel = new JLabel(" ");
        this.statusLabel.setFont(new Font("SansSerif", 2, 11));
        this.statusLabel.setForeground(Color.BLUE);
        jPanel4.add(this.statusLabel);
        jPanel3.add((Component)jPanel4, "West");
        JPanel jPanel5 = new JPanel(new FlowLayout(2));
        jPanel5.add(new JLabel("Var Min \u0394V:"));
        this.varDeltaVSpinner = new JSpinner(new SpinnerNumberModel(0.06, 0.0, 5.0, 0.01));
        this.varDeltaVSpinner.setPreferredSize(new Dimension(70, 25));
        this.varDeltaVSpinner.addChangeListener(changeEvent -> {
            this.minVarDeltaV = (Double)this.varDeltaVSpinner.getValue();
            this.imagePanel.repaint();
        });
        jPanel5.add(this.varDeltaVSpinner);
        jPanel5.add(Box.createHorizontalStrut(10));
        jPanel5.add(new JLabel("FOV (arcmin):"));
        this.fovSpinner = new JSpinner(new SpinnerNumberModel(60.0, 1.0, 120.0, 1.0));
        this.fovSpinner.setPreferredSize(new Dimension(70, 25));
        this.fovSpinner.addChangeListener(changeEvent -> {
            this.currentFOV = (Double)this.fovSpinner.getValue();
        });
        jPanel5.add(this.fovSpinner);
        jPanel5.add(Box.createHorizontalStrut(10));
        jPanel5.add(new JLabel("Mag Limit:"));
        this.magLimitSpinner = new JSpinner(new SpinnerNumberModel(14.5, 8.0, 20.0, 0.5));
        this.magLimitSpinner.setPreferredSize(new Dimension(70, 25));
        this.magLimitSpinner.addChangeListener(changeEvent -> {
            this.currentMagLimit = (Double)this.magLimitSpinner.getValue();
        });
        jPanel5.add(this.magLimitSpinner);
        jPanel3.add((Component)jPanel5, "East");
        jPanel2.add((Component)jPanel3, "North");
        JPanel jPanel6 = new JPanel(new FlowLayout(0, 10, 2));
        jPanel6.setBorder(BorderFactory.createEmptyBorder(0, 5, 0, 5));
        this.mouseInfoLabel = new JLabel("x: -, y: -, intensity: -");
        this.mouseInfoLabel.setFont(new Font("Monospaced", 0, 12));
        this.mouseInfoLabel.setForeground(new Color(0, 150, 0));
        jPanel6.add(this.mouseInfoLabel);
        jPanel.add((Component)jPanel6, "North");
        this.detailsArea = new JTextArea(3, 120);
        this.detailsArea.setEditable(false);
        this.detailsArea.setFont(new Font("Monospaced", 0, 11));
        this.detailsArea.setBorder(BorderFactory.createTitledBorder("Star Details"));
        this.detailsArea.setLineWrap(true);
        this.detailsArea.setWrapStyleWord(true);
        JScrollPane jScrollPane = new JScrollPane(this.detailsArea);
        jScrollPane.setVerticalScrollBarPolicy(20);
        jScrollPane.setPreferredSize(new Dimension(0, 85));
        jPanel2.add((Component)jScrollPane, "Center");
        jPanel.add((Component)jPanel2, "South");
        return jPanel;
    }

    private void refreshTargetList() {
        this.tableModel.setRowCount(0);
        List<TargetNightEntry> list = this.scanTargetNights();
        list.sort((targetNightEntry, targetNightEntry2) -> {
            int n = this.getStatusPriority(targetNightEntry.status) - this.getStatusPriority(targetNightEntry2.status);
            if (n != 0) {
                return n;
            }
            return targetNightEntry2.night.compareTo(targetNightEntry.night);
        });
        for (TargetNightEntry targetNightEntry3 : list) {
            this.tableModel.addRow(new Object[]{targetNightEntry3.target, targetNightEntry3.night.format(DateTimeFormatter.ofPattern("yyyy-MM-dd")), targetNightEntry3.status, targetNightEntry3.imageCount});
        }
    }

    private int getStatusPriority(String string) {
        switch (string) {
            case "Not Done": {
                return 0;
            }
            case "Done": {
                return 1;
            }
            case "Bad Obs": {
                return 2;
            }
        }
        return 3;
    }

    private List<TargetNightEntry> scanTargetNights() {
        ArrayList<TargetNightEntry> arrayList = new ArrayList<TargetNightEntry>();
        File file2 = new File(this.workspaceRoot, "data");
        if (!file2.exists() || !file2.isDirectory()) {
            return arrayList;
        }
        Map<String, String> map = this.loadPhotometryStatus();
        for (File file3 : file2.listFiles()) {
            if (!file3.isDirectory()) continue;
            try {
                File[] fileArray;
                LocalDate localDate = LocalDate.parse(file3.getName(), DateTimeFormatter.ofPattern("yyyyMMdd"));
                File file4 = new File(file3, "Light");
                if (!file4.exists() || !file4.isDirectory() || (fileArray = file4.listFiles(File::isDirectory)) == null) continue;
                for (File file5 : fileArray) {
                    int n;
                    String string2 = file5.getName();
                    File file6 = new File(file5, "Reduced_images");
                    if (!file6.exists() || !file6.isDirectory()) continue;
                    File[] fileArray2 = file6.listFiles((file, string) -> string.toLowerCase().endsWith("_cal.fits"));
                    int n2 = n = fileArray2 != null ? fileArray2.length : 0;
                    if (n == 0) continue;
                    String string3 = file3.getName() + "_" + string2;
                    String string4 = map.getOrDefault(string3, STATUS_NOT_DONE);
                    arrayList.add(new TargetNightEntry(string2, localDate, string4, n));
                }
            }
            catch (Exception exception) {
                // empty catch block
            }
        }
        return arrayList;
    }

    private Map<String, String> loadPhotometryStatus() {
        HashMap<String, String> hashMap = new HashMap<String, String>();
        File file = new File(this.workspaceRoot, STATUS_FILE);
        if (!file.exists()) {
            return hashMap;
        }
        try (BufferedReader bufferedReader = new BufferedReader(new FileReader(file));){
            String string;
            while ((string = bufferedReader.readLine()) != null) {
                String[] stringArray;
                if ((string = string.trim()).isEmpty() || string.startsWith("#") || (stringArray = string.split("\\t")).length != 2) continue;
                hashMap.put(stringArray[0], stringArray[1]);
            }
        }
        catch (Exception exception) {
            System.err.println("Error loading photometry status: " + exception.getMessage());
        }
        return hashMap;
    }

    private void savePhotometryStatus(String string, String string2, String string3) {
        Map<String, String> map = this.loadPhotometryStatus();
        String string4 = string + "_" + string2;
        map.put(string4, string3);
        File file = new File(this.workspaceRoot, STATUS_FILE);
        try (PrintWriter printWriter = new PrintWriter(new FileWriter(file));){
            printWriter.println("# Photometry Processing Status");
            printWriter.println("# Format: nightDirName_targetName<TAB>status");
            for (Map.Entry<String, String> entry : map.entrySet()) {
                printWriter.println(entry.getKey() + "\t" + entry.getValue());
            }
        }
        catch (Exception exception) {
            System.err.println("Error saving photometry status: " + exception.getMessage());
        }
    }

    private void onTargetSelected() {
        File[] fileArray;
        File file2;
        int n = this.targetTable.getSelectedRow();
        if (n < 0) {
            this.imagePanel.setImage(null, null, null);
            this.imagePanel.setMessage("Select a target-night to view image");
            this.currentFitsFile = null;
            this.currentFitsHeader = null;
            this.currentImage = null;
            this.referenceStars = null;
            return;
        }
        String string2 = (String)this.tableModel.getValueAt(n, 0);
        String string3 = (String)this.tableModel.getValueAt(n, 1);
        LocalDate localDate = LocalDate.parse(string3, DateTimeFormatter.ofPattern("yyyy-MM-dd"));
        String string4 = localDate.format(DateTimeFormatter.ofPattern("yyyyMMdd"));
        this.currentTargetDir = file2 = new File(this.workspaceRoot, "data/" + string4 + "/Light/" + string2);
        File file3 = null;
        File file4 = new File(file2, "Reduced_images");
        if (file4.exists() && (fileArray = file4.listFiles((file, string) -> string.toLowerCase().endsWith("_cal.fits"))) != null && fileArray.length > 0) {
            file3 = fileArray[0];
        }
        if (file3 == null && (fileArray = file2.listFiles((file, string) -> string.toLowerCase().endsWith(".fits"))) != null && fileArray.length > 0) {
            file3 = fileArray[0];
        }
        if (file3 != null) {
            this.loadAndDisplayImage(file3);
        } else {
            this.imagePanel.setImage(null, null, null);
            this.imagePanel.setMessage("No images found for " + string2 + " on " + string3);
            this.currentFitsFile = null;
            this.currentFitsHeader = null;
            this.currentImage = null;
            this.referenceStars = null;
        }
    }

    private void loadAndDisplayImage(File file) {
        try {
            this.imagePanel.setMessage("Loading " + file.getName() + "...");
            Fits fits = new Fits(file);
            ImageHDU imageHDU = (ImageHDU)fits.getHDU(0);
            Object object = ((ImageData)imageHDU.getData()).getData();
            Header header = imageHDU.getHeader();
            int[][] nArray = this.convertToIntArray(object);
            int n = Integer.MAX_VALUE;
            int n2 = Integer.MIN_VALUE;
            for (int i = 0; i < nArray.length; ++i) {
                for (int j = 0; j < nArray[i].length; ++j) {
                    n = Math.min(n, nArray[i][j]);
                    n2 = Math.max(n2, nArray[i][j]);
                }
            }
            System.out.println("IMAGE LOADED: " + file.getName() + " - dimensions=" + nArray.length + "x" + nArray[0].length + ", min=" + n + ", max=" + n2);
            BufferedImage bufferedImage = this.createImageFromData(nArray);
            this.currentFitsFile = file;
            this.currentFitsHeader = header;
            this.currentImage = bufferedImage;
            this.currentImageData = nArray;
            this.referenceStars = null;
            this.calculateFieldOfView(header, nArray);
            this.imagePanel.setImage(bufferedImage, header, this.referenceStars);
            this.imagePanel.setMessage(null);
            fits.close();
            SwingUtilities.invokeLater(() -> {
                this.onGetVSD();
                this.generateRadialProfile();
            });
        }
        catch (Exception exception) {
            this.imagePanel.setImage(null, null, null);
            this.imagePanel.setMessage("Error loading image: " + exception.getMessage());
            this.currentFitsFile = null;
            this.currentFitsHeader = null;
            this.currentImage = null;
            this.referenceStars = null;
            exception.printStackTrace();
        }
    }

    private int[][] convertToIntArray(Object object) {
        if (object instanceof short[][]) {
            short[][] sArray = (short[][])object;
            int[][] nArray = new int[sArray.length][sArray[0].length];
            int n = Integer.MAX_VALUE;
            int n2 = Integer.MIN_VALUE;
            for (int i = 0; i < sArray.length; ++i) {
                for (int j = 0; j < sArray[i].length; ++j) {
                    int n3;
                    nArray[i][j] = n3 = sArray[i][j] & 0xFFFF;
                    if (n3 < n) {
                        n = n3;
                    }
                    if (n3 <= n2) continue;
                    n2 = n3;
                }
            }
            System.out.println("FITS DATA: Converted short[][] to int[][] - min=" + n + ", max=" + n2 + ", dimensions=" + sArray.length + "x" + sArray[0].length);
            return nArray;
        }
        if (object instanceof int[][]) {
            int[][] nArray = (int[][])object;
            int n = Integer.MAX_VALUE;
            int n4 = Integer.MIN_VALUE;
            for (int i = 0; i < nArray.length; ++i) {
                for (int j = 0; j < nArray[i].length; ++j) {
                    if (nArray[i][j] < n) {
                        n = nArray[i][j];
                    }
                    if (nArray[i][j] <= n4) continue;
                    n4 = nArray[i][j];
                }
            }
            System.out.println("FITS DATA: int[][] data - min=" + n + ", max=" + n4 + ", dimensions=" + nArray.length + "x" + nArray[0].length);
            return nArray;
        }
        if (object instanceof float[][]) {
            float[][] fArray = (float[][])object;
            int[][] nArray = new int[fArray.length][fArray[0].length];
            float f = Float.MAX_VALUE;
            float f2 = Float.MIN_VALUE;
            for (int i = 0; i < fArray.length; ++i) {
                for (int j = 0; j < fArray[i].length; ++j) {
                    nArray[i][j] = (int)fArray[i][j];
                    if (fArray[i][j] < f) {
                        f = fArray[i][j];
                    }
                    if (!(fArray[i][j] > f2)) continue;
                    f2 = fArray[i][j];
                }
            }
            System.out.println("FITS DATA: Converted float[][] to int[][] - min=" + f + ", max=" + f2 + ", dimensions=" + fArray.length + "x" + fArray[0].length);
            return nArray;
        }
        System.err.println("FITS DATA ERROR: Unsupported data type: " + String.valueOf(object.getClass()));
        throw new IllegalArgumentException("Unsupported data type: " + String.valueOf(object.getClass()));
    }

    private BufferedImage createImageFromData(int[][] nArray) {
        int n;
        int n2;
        int n3 = nArray.length;
        int n4 = nArray[0].length;
        int[] nArray2 = new int[n4 * n3];
        int n5 = 0;
        for (n2 = 0; n2 < n3; ++n2) {
            for (n = 0; n < n4; ++n) {
                nArray2[n5++] = nArray[n2][n];
            }
        }
        Arrays.sort(nArray2);
        n2 = (int)((double)nArray2.length * 0.01);
        n = (int)((double)nArray2.length * 0.99);
        int n6 = nArray2[n2];
        int n7 = nArray2[n];
        if (n7 <= n6) {
            n7 = n6 + 1;
        }
        BufferedImage bufferedImage = new BufferedImage(n4, n3, 10);
        for (int i = 0; i < n3; ++i) {
            for (int j = 0; j < n4; ++j) {
                int n8 = nArray[i][j];
                int n9 = n8 <= n6 ? 0 : (n8 >= n7 ? 255 : (int)(255.0 * (double)(n8 - n6) / (double)(n7 - n6)));
                n9 = Math.max(0, Math.min(255, n9));
                int n10 = n9 << 16 | n9 << 8 | n9;
                bufferedImage.setRGB(j, i, n10);
            }
        }
        return bufferedImage;
    }

    private void onGetVSD() {
        if (this.currentFitsFile == null || this.currentFitsHeader == null) {
            JOptionPane.showMessageDialog(this, "No image loaded", "Error", 0);
            return;
        }
        int n = this.targetTable.getSelectedRow();
        if (n < 0) {
            return;
        }
        final String string = (String)this.tableModel.getValueAt(n, 0);
        final double d = this.currentFOV;
        final double d2 = this.currentMagLimit;
        System.out.println("=== VSX/VSD Query Starting ===");
        System.out.println("Target: " + string);
        System.out.println("FOV: " + d + " arcmin");
        System.out.println("Mag Limit: " + d2);
        this.statusLabel.setText("Searching database...");
        this.statusLabel.setForeground(Color.BLUE);
        SwingWorker<Void, Void> swingWorker = new SwingWorker<Void, Void>(this){
            final /* synthetic */ TargetPhotometryPanel this$0;
            {
                this.this$0 = targetPhotometryPanel;
            }

            @Override
            protected Void doInBackground() {
                try {
                    this.this$0.targetInfo = this.this$0.fetchVSXTargetInfo(string);
                    this.this$0.variableStarsInField = new ArrayList<TargetInfo>();
                    if (this.this$0.targetInfo != null) {
                        System.out.println("VSX Target Info:");
                        System.out.println("  Name: " + this.this$0.targetInfo.name);
                        System.out.println("  AUID: " + this.this$0.targetInfo.auid);
                        System.out.println("  RA: " + this.this$0.targetInfo.ra + "\u00b0");
                        System.out.println("  Dec: " + this.this$0.targetInfo.dec + "\u00b0");
                        System.out.println("  Type: " + this.this$0.targetInfo.variabilityType);
                        System.out.println("  Period: " + this.this$0.targetInfo.period + " days");
                        System.out.println("  Mag range: " + this.this$0.targetInfo.maxMag + " to " + this.this$0.targetInfo.minMag);
                        System.out.println("  Constellation: " + this.this$0.targetInfo.constellation);
                        this.this$0.variableStarsInField.add(this.this$0.targetInfo);
                    }
                    this.this$0.referenceStars = this.this$0.fetchAAVSOVSD(string, null, null, d, d2);
                    Double d3 = this.this$0.getCenterRA(this.this$0.currentFitsHeader);
                    Double d22 = this.this$0.getCenterDec(this.this$0.currentFitsHeader);
                    if (d3 != null && d22 != null) {
                        System.out.println("Field center: RA=" + d3 + ", Dec=" + d22);
                        this.this$0.fetchVSXFieldStars(d3, d22, d);
                    }
                    if (this.this$0.referenceStars == null || this.this$0.referenceStars.isEmpty()) {
                        if (d3 != null && d22 != null) {
                            System.out.println("Target name not found, trying coordinates: RA=" + d3 + ", Dec=" + d22);
                            this.this$0.referenceStars = this.this$0.fetchAAVSOVSD(null, d3, d22, d, d2);
                        } else {
                            SwingUtilities.invokeLater(() -> JOptionPane.showMessageDialog(this.this$0, "Target name not recognized and no WCS coordinates found in image", "VSD Error", 2));
                            return null;
                        }
                    }
                    if (this.this$0.referenceStars == null || this.this$0.referenceStars.isEmpty()) {
                        SwingUtilities.invokeLater(() -> JOptionPane.showMessageDialog(this.this$0, "No reference stars found or error fetching VSD data", "VSD Error", 2));
                        return null;
                    }
                }
                catch (Exception exception) {
                    exception.printStackTrace();
                    SwingUtilities.invokeLater(() -> JOptionPane.showMessageDialog(this.this$0, "Error during database query: " + exception.getMessage(), "Query Error", 0));
                }
                return null;
            }

            @Override
            protected void done() {
                if (this.this$0.referenceStars != null && !this.this$0.referenceStars.isEmpty()) {
                    this.this$0.calculateMagnitudeLimit();
                    this.this$0.imagePanel.setImage(this.this$0.currentImage, this.this$0.currentFitsHeader, this.this$0.referenceStars);
                    String string2 = "Loaded " + this.this$0.referenceStars.size() + " reference stars";
                    if (this.this$0.variableStarsInField != null && !this.this$0.variableStarsInField.isEmpty()) {
                        string2 = string2 + "\nFound " + this.this$0.variableStarsInField.size() + " variable star(s) in field";
                    }
                    if (this.this$0.targetInfo != null) {
                        string2 = string2 + "\n\nTarget: " + this.this$0.targetInfo.name + " (" + this.this$0.targetInfo.variabilityType + ")\nMag range: " + this.this$0.targetInfo.maxMag + " to " + this.this$0.targetInfo.minMag;
                        if (this.this$0.targetInfo.period > 0.0) {
                            string2 = string2 + "\nPeriod: " + String.format("%.5f days", this.this$0.targetInfo.period);
                        }
                    }
                    this.this$0.statusLabel.setText("Query complete");
                    this.this$0.statusLabel.setForeground(new Color(0, 128, 0));
                    System.out.println("=== VSX/VSD Query Complete ===");
                } else {
                    this.this$0.statusLabel.setText("Query failed");
                    this.this$0.statusLabel.setForeground(Color.RED);
                }
            }
        };
        swingWorker.execute();
    }

    private void calculateFieldOfView(Header header, int[][] nArray) {
        try {
            int n = nArray.length;
            int n2 = nArray[0].length;
            Double d = header.getDoubleValue("CD1_1");
            Double d2 = header.getDoubleValue("CD1_2");
            Double d3 = header.getDoubleValue("CD2_1");
            Double d4 = header.getDoubleValue("CD2_2");
            double d5 = 0.0;
            if (d != null && d4 != null) {
                d5 = d3 != null ? Math.sqrt(d * d + d3 * d3) : Math.abs(d);
            } else {
                Double d6 = header.getDoubleValue("CDELT1");
                Double d7 = header.getDoubleValue("CDELT2");
                if (d6 != null && d7 != null) {
                    d5 = Math.sqrt(d6 * d6 + d7 * d7);
                }
            }
            if (d5 > 0.0) {
                double d8 = (double)n * d5 * 60.0;
                double d9 = (double)n2 * d5 * 60.0;
                this.currentFOV = Math.max(d8, d9);
                System.out.println("Calculated FOV: " + String.format("%.2f", this.currentFOV) + " arcmin (X=" + String.format("%.2f", d8) + ", Y=" + String.format("%.2f", d9) + ")");
                SwingUtilities.invokeLater(() -> this.fovSpinner.setValue(this.currentFOV));
            } else {
                System.out.println("Could not calculate FOV - no WCS scale information");
            }
        }
        catch (Exception exception) {
            System.err.println("Error calculating FOV: " + exception.getMessage());
        }
    }

    private void calculateMagnitudeLimit() {
        if (this.referenceStars == null || this.referenceStars.isEmpty()) {
            return;
        }
        try {
            double d = Double.MIN_VALUE;
            for (ReferenceStar referenceStar : this.referenceStars) {
                if (!(referenceStar.mag > d)) continue;
                d = referenceStar.mag;
            }
            if (d > Double.MIN_VALUE) {
                this.currentMagLimit = d;
                System.out.println("Calculated Mag Limit: " + String.format("%.2f", this.currentMagLimit) + " (from " + this.referenceStars.size() + " comp stars)");
                SwingUtilities.invokeLater(() -> this.magLimitSpinner.setValue(this.currentMagLimit));
            }
        }
        catch (Exception exception) {
            System.err.println("Error calculating mag limit: " + exception.getMessage());
        }
    }

    private JPanel createRadialProfileTab() {
        JPanel jPanel = new JPanel(new BorderLayout(10, 10));
        jPanel.setBorder(BorderFactory.createEmptyBorder(10, 10, 10, 10));
        JLabel jLabel = new JLabel("Load an image to generate radial profile", 0);
        jLabel.setFont(new Font("SansSerif", 2, 12));
        jPanel.add((Component)jLabel, "Center");
        return jPanel;
    }

    private JPanel createPhotometryResultsTab() {
        this.photometryResultsPanel = new JPanel(new BorderLayout(10, 10));
        this.photometryResultsPanel.setBorder(BorderFactory.createEmptyBorder(10, 10, 10, 10));
        JLabel jLabel = new JLabel("Click 'Do Phot' to perform photometry", 0);
        jLabel.setFont(new Font("SansSerif", 2, 12));
        this.photometryResultsPanel.add((Component)jLabel, "Center");
        return this.photometryResultsPanel;
    }

    private JPanel createLightcurvesTab() {
        this.lightcurvePanel = new JPanel(new BorderLayout(10, 10));
        this.lightcurvePanel.setBorder(BorderFactory.createEmptyBorder(10, 10, 10, 10));
        JPanel jPanel = new JPanel(new FlowLayout(0));
        jPanel.add(new JLabel("Select Variable:"));
        this.variableSelector = new JComboBox();
        this.variableSelector.setPreferredSize(new Dimension(200, 25));
        this.variableSelector.addActionListener(actionEvent -> this.updateLightcurve());
        jPanel.add(this.variableSelector);
        this.lightcurvePanel.add((Component)jPanel, "North");
        JLabel jLabel = new JLabel("Run batch photometry to generate lightcurves", 0);
        jLabel.setFont(new Font("SansSerif", 2, 12));
        this.lightcurvePanel.add((Component)jLabel, "Center");
        return this.lightcurvePanel;
    }

    private void updateLightcurve() {
        Component[] componentArray;
        if (this.photometryResults == null || this.photometryResults.isEmpty() || this.variableSelector.getSelectedItem() == null) {
            return;
        }
        String string = (String)this.variableSelector.getSelectedItem();
        ArrayList<LightcurvePoint> arrayList = new ArrayList<LightcurvePoint>();
        ArrayList<LightcurvePoint> arrayList2 = new ArrayList<LightcurvePoint>();
        for (PhotometryResult componentArray2 : this.photometryResults) {
            double d;
            double d2;
            if (Double.isNaN(componentArray2.calibratedMag) || !((d2 = this.getHJDFromImageId(componentArray2.imageId)) > 0.0)) continue;
            if (componentArray2.name.equals(string) && componentArray2.type.equals("Variable")) {
                d = Math.sqrt(componentArray2.magError * componentArray2.magError + componentArray2.zpUncertainty * componentArray2.zpUncertainty);
                arrayList.add(new LightcurvePoint(d2, componentArray2.calibratedMag, d));
                continue;
            }
            if (!componentArray2.type.equals("Comparison") || this.checkStar == null || !componentArray2.name.equals(this.checkStar.label) && !componentArray2.name.equals(this.checkStar.name)) continue;
            d = Math.sqrt(componentArray2.magError * componentArray2.magError + componentArray2.zpUncertainty * componentArray2.zpUncertainty);
            arrayList2.add(new LightcurvePoint(d2, componentArray2.calibratedMag, d));
        }
        if (arrayList.isEmpty()) {
            return;
        }
        LightcurvePlotPanel lightcurvePlotPanel = new LightcurvePlotPanel(string, arrayList, arrayList2);
        for (Component component : componentArray = this.lightcurvePanel.getComponents()) {
            if (component == this.lightcurvePanel.getComponent(0)) continue;
            this.lightcurvePanel.remove(component);
        }
        this.lightcurvePanel.add((Component)lightcurvePlotPanel, "Center");
        this.lightcurvePanel.revalidate();
        this.lightcurvePanel.repaint();
    }

    private double getHJDFromImageId(String string) {
        try {
            File file;
            File file2 = new File(this.currentTargetDir, string);
            if (!file2.exists() && (file = new File(this.currentTargetDir, "Reduced_images")).exists()) {
                file2 = new File(file, string);
            }
            if (file2.exists()) {
                file = new Fits(file2);
                ImageHDU imageHDU = (ImageHDU)file.readHDU();
                Header header = imageHDU.getHeader();
                String string2 = header.getStringValue("HJD_UTC");
                if (string2 != null && !string2.trim().isEmpty()) {
                    double d = Double.parseDouble(string2);
                    file.close();
                    return d;
                }
                string2 = header.getStringValue("HJD");
                if (string2 != null && !string2.trim().isEmpty()) {
                    double d = Double.parseDouble(string2);
                    file.close();
                    return d;
                }
                String string3 = header.getStringValue("BJD_TDB");
                if (string3 != null && !string3.trim().isEmpty()) {
                    double d = Double.parseDouble(string3);
                    file.close();
                    return d;
                }
                String string4 = header.getStringValue("JD_UTC");
                if (string4 != null && !string4.trim().isEmpty()) {
                    double d = Double.parseDouble(string4);
                    file.close();
                    return d;
                }
                string4 = header.getStringValue("JD");
                if (string4 != null && !string4.trim().isEmpty()) {
                    double d = Double.parseDouble(string4);
                    file.close();
                    return d;
                }
                String string5 = header.getStringValue("DATE-AVG");
                if (string5 == null || string5.trim().isEmpty()) {
                    string5 = header.getStringValue("DATE-OBS");
                }
                if (string5 != null && !string5.trim().isEmpty()) {
                    double d = this.isoToJD(string5);
                    file.close();
                    return d;
                }
                file.close();
            }
        }
        catch (Exception exception) {
            System.err.println("Error extracting HJD from " + string + ": " + exception.getMessage());
        }
        return 0.0;
    }

    private double isoToJD(String string) {
        try {
            String[] stringArray = string.split("T");
            String[] stringArray2 = stringArray[0].split("-");
            String[] stringArray3 = stringArray[1].split(":");
            int n = Integer.parseInt(stringArray2[0]);
            int n2 = Integer.parseInt(stringArray2[1]);
            int n3 = Integer.parseInt(stringArray2[2]);
            double d = Double.parseDouble(stringArray3[0]);
            double d2 = Double.parseDouble(stringArray3[1]);
            double d3 = Double.parseDouble(stringArray3[2]);
            int n4 = (14 - n2) / 12;
            int n5 = n + 4800 - n4;
            int n6 = n2 + 12 * n4 - 3;
            double d4 = n3 + (153 * n6 + 2) / 5 + 365 * n5 + n5 / 4 - n5 / 100 + n5 / 400 - 32045;
            double d5 = (d - 12.0) / 24.0 + d2 / 1440.0 + d3 / 86400.0;
            return d4 + d5;
        }
        catch (Exception exception) {
            System.err.println("Error parsing date: " + string + ": " + exception.getMessage());
            return 0.0;
        }
    }

    private void generateRadialProfile() {
        RadialProfile radialProfile;
        if (this.currentImageData == null || this.currentFitsHeader == null) {
            return;
        }
        StarLocation starLocation = this.findSuitableStar(this.currentImageData, 1000, 30000);
        if (starLocation == null) {
            this.radialProfileStar = null;
            this.imagePanel.repaint();
            this.radialProfilePanel.removeAll();
            JLabel jLabel = new JLabel("No suitable star found for radial profile", 0);
            jLabel.setFont(new Font("SansSerif", 2, 12));
            this.radialProfilePanel.add((Component)jLabel, "Center");
            this.radialProfilePanel.revalidate();
            this.radialProfilePanel.repaint();
            return;
        }
        System.out.println("Found suitable star at (" + starLocation.x + ", " + starLocation.y + ") with peak " + starLocation.peakValue);
        this.radialProfileStar = starLocation;
        this.imagePanel.repaint();
        this.currentRadialProfile = radialProfile = this.calculateRadialProfile(this.currentImageData, starLocation.x, starLocation.y, 50);
        this.updateRadialProfilePanel(radialProfile, starLocation);
    }

    private void updateRadialProfilePanel(RadialProfile radialProfile, StarLocation starLocation) {
        this.radialProfilePanel.removeAll();
        this.radialProfilePanel.setLayout(new BorderLayout(10, 10));
        RadialProfilePlot radialProfilePlot = new RadialProfilePlot(radialProfile);
        radialProfilePlot.setPreferredSize(new Dimension(700, 500));
        this.radialProfilePanel.add((Component)radialProfilePlot, "Center");
        JPanel jPanel = new JPanel(new FlowLayout(1, 20, 10));
        jPanel.setBorder(BorderFactory.createTitledBorder("Aperture Settings"));
        jPanel.add(new JLabel("Star Location:"));
        jPanel.add(new JLabel(String.format("(%d, %d)", starLocation.x, starLocation.y)));
        jPanel.add(Box.createHorizontalStrut(20));
        jPanel.add(new JLabel("Peak Value:"));
        jPanel.add(new JLabel(String.format("%d counts", starLocation.peakValue)));
        jPanel.add(Box.createHorizontalStrut(20));
        jPanel.add(new JLabel("Object Radius:"));
        JSpinner jSpinner = new JSpinner(new SpinnerNumberModel(radialProfile.objectRadius, 1, 50, 1));
        jPanel.add(jSpinner);
        jPanel.add(Box.createHorizontalStrut(10));
        jPanel.add(new JLabel("Sky Inner Radius:"));
        JSpinner jSpinner2 = new JSpinner(new SpinnerNumberModel(radialProfile.skyInnerRadius, 1, 100, 1));
        jPanel.add(jSpinner2);
        jPanel.add(Box.createHorizontalStrut(10));
        jPanel.add(new JLabel("Sky Outer Radius:"));
        JSpinner jSpinner3 = new JSpinner(new SpinnerNumberModel(radialProfile.skyOuterRadius, 1, 100, 1));
        jPanel.add(jSpinner3);
        jPanel.add(Box.createHorizontalStrut(20));
        JButton jButton = new JButton("Update Plot");
        jButton.addActionListener(actionEvent -> {
            radialProfile.objectRadius = (Integer)jSpinner.getValue();
            radialProfile.skyInnerRadius = (Integer)jSpinner2.getValue();
            radialProfile.skyOuterRadius = (Integer)jSpinner3.getValue();
            this.currentRadialProfile = radialProfile;
            radialProfilePlot.updateProfile(radialProfile);
        });
        jPanel.add(jButton);
        this.radialProfilePanel.add((Component)jPanel, "South");
        this.radialProfilePanel.revalidate();
        this.radialProfilePanel.repaint();
    }

    private void onRadialProfile() {
        if (this.currentImageData == null || this.currentFitsHeader == null) {
            JOptionPane.showMessageDialog(this, "No image loaded", "Error", 0);
            return;
        }
        StarLocation starLocation = this.findSuitableStar(this.currentImageData, 1000, 30000);
        if (starLocation == null) {
            JOptionPane.showMessageDialog(this, "No suitable star found (looking for 1000-30000 counts)", "No Star Found", 2);
            return;
        }
        System.out.println("Found suitable star at (" + starLocation.x + ", " + starLocation.y + ") with peak " + starLocation.peakValue);
        RadialProfile radialProfile = this.calculateRadialProfile(this.currentImageData, starLocation.x, starLocation.y, 50);
        this.showRadialProfileDialog(radialProfile, starLocation);
    }

    private StarLocation findSuitableStar(int[][] nArray, int n, int n2) {
        int n3 = nArray.length;
        int n4 = nArray[0].length;
        ArrayList<StarLocation> arrayList = new ArrayList<StarLocation>();
        for (int i = 20; i < n3 - 20; ++i) {
            for (int j = 20; j < n4 - 20; ++j) {
                int n5 = nArray[i][j];
                if (n5 < n || n5 > n2) continue;
                boolean bl = true;
                for (int k = -2; k <= 2 && bl; ++k) {
                    for (int i2 = -2; i2 <= 2 && bl; ++i2) {
                        if (i2 == 0 && k == 0 || nArray[i + k][j + i2] <= n5) continue;
                        bl = false;
                    }
                }
                if (!bl) continue;
                arrayList.add(new StarLocation(j, i, n5));
            }
        }
        if (arrayList.isEmpty()) {
            return null;
        }
        arrayList.sort((starLocation, starLocation2) -> Integer.compare(starLocation2.peakValue, starLocation.peakValue));
        return (StarLocation)arrayList.get(0);
    }

    private RadialProfile calculateRadialProfile(int[][] nArray, int n, int n2, int n3) {
        int n4;
        int n5 = nArray.length;
        int n6 = nArray[0].length;
        double[] dArray = new double[n3 + 1];
        int[] nArray2 = new int[n3 + 1];
        for (int i = Math.max(0, n2 - n3); i <= Math.min(n5 - 1, n2 + n3); ++i) {
            for (n4 = Math.max(0, n - n3); n4 <= Math.min(n6 - 1, n + n3); ++n4) {
                double d = n4 - n;
                double d2 = i - n2;
                double d3 = Math.sqrt(d * d + d2 * d2);
                int n7 = (int)Math.round(d3);
                if (n7 > n3) continue;
                int n8 = n7;
                dArray[n8] = dArray[n8] + (double)nArray[i][n4];
                int n9 = n7;
                nArray2[n9] = nArray2[n9] + 1;
            }
        }
        double[] dArray2 = new double[n3 + 1];
        for (n4 = 0; n4 <= n3; ++n4) {
            if (nArray2[n4] <= 0) continue;
            dArray2[n4] = dArray[n4] / (double)nArray2[n4];
        }
        // Use FWHM-based aperture sizing like AstroImageJ
        double fwhm = this.calculateFWHM(dArray2);
        n4 = (int)Math.round(fwhm * 0.85);  // Object aperture = 0.85 * FWHM
        int n10 = (int)Math.round(fwhm * 1.7);  // Sky inner = 1.7 * FWHM
        int n11 = (int)Math.round(fwhm * 3.4);  // Sky outer = 3.4 * FWHM
        // Ensure minimum reasonable values
        n4 = Math.max(3, n4);
        n10 = Math.max(n4 + 2, n10);
        n11 = Math.max(n10 + 2, n11);
        System.out.println(String.format("FWHM-based apertures: FWHM=%.2f, object=%d, sky_inner=%d, sky_outer=%d", 
            fwhm, n4, n10, n11));
        return new RadialProfile(dArray2, n4, n10, n11, n, n2);
    }

    private int suggestObjectRadius(double[] dArray) {
        if (dArray.length < 2) {
            return 5;
        }
        double d = dArray[0];
        double d2 = d * 0.1;
        for (int i = 1; i < dArray.length; ++i) {
            if (!(dArray[i] < d2)) continue;
            return Math.max(3, i + 2);
        }
        return Math.min(10, dArray.length / 2);
    }

    /**
     * Calculate Full Width at Half Maximum (FWHM) from radial profile.
     * Uses the same approach as AstroImageJ Seeing_Profile.
     * 
     * @param profile Array of average intensities at each radius
     * @return FWHM in pixels
     */
    private double calculateFWHM(double[] profile) {
        if (profile.length < 2) {
            return 10.0; // Default fallback
        }
        
        // Find peak intensity (at center)
        double peak = profile[0];
        
        // Estimate background from the tail of the profile
        double background = 0.0;
        int nBackgroundPoints = 0;
        int startIdx = Math.max(1, profile.length * 2 / 3);
        for (int i = startIdx; i < profile.length; i++) {
            if (profile[i] > 0) {
                background += profile[i];
                nBackgroundPoints++;
            }
        }
        if (nBackgroundPoints > 0) {
            background /= nBackgroundPoints;
        }
        
        // Calculate half-maximum level
        double halfMax = background + (peak - background) / 2.0;
        
        // Find radius where intensity drops to half-maximum
        double halfMaxRadius = 0.0;
        for (int i = 1; i < profile.length; i++) {
            if (profile[i] <= halfMax) {
                // Interpolate between this point and previous
                if (i > 0 && profile[i-1] > halfMax) {
                    double fraction = (halfMax - profile[i]) / (profile[i-1] - profile[i]);
                    halfMaxRadius = i - fraction;
                } else {
                    halfMaxRadius = i;
                }
                break;
            }
        }
        
        // FWHM = 2 * radius at half maximum
        double fwhm = 2.0 * halfMaxRadius;
        
        // Sanity check: FWHM should be reasonable (between 2 and 50 pixels for typical seeing)
        if (fwhm < 2.0) {
            fwhm = 2.0;
        } else if (fwhm > 50.0) {
            fwhm = 50.0;
        }
        
        return fwhm;
    }

    private void showRadialProfileDialog(RadialProfile radialProfile, StarLocation starLocation) {
        JDialog jDialog = new JDialog((Frame)SwingUtilities.getWindowAncestor(this), "Radial Profile", true);
        jDialog.setLayout(new BorderLayout(10, 10));
        RadialProfilePlot radialProfilePlot = new RadialProfilePlot(radialProfile);
        radialProfilePlot.setPreferredSize(new Dimension(600, 400));
        jDialog.add((Component)radialProfilePlot, "Center");
        JPanel jPanel = new JPanel(new GridLayout(6, 2, 5, 5));
        jPanel.setBorder(BorderFactory.createTitledBorder("Aperture Settings"));
        jPanel.add(new JLabel("Star Location:"));
        jPanel.add(new JLabel(String.format("(%d, %d)", starLocation.x, starLocation.y)));
        jPanel.add(new JLabel("Peak Value:"));
        jPanel.add(new JLabel(String.format("%d counts", starLocation.peakValue)));
        jPanel.add(new JLabel("Object Radius:"));
        JSpinner jSpinner = new JSpinner(new SpinnerNumberModel(radialProfile.objectRadius, 1, 50, 1));
        jPanel.add(jSpinner);
        jPanel.add(new JLabel("Sky Inner Radius:"));
        JSpinner jSpinner2 = new JSpinner(new SpinnerNumberModel(radialProfile.skyInnerRadius, 1, 100, 1));
        jPanel.add(jSpinner2);
        jPanel.add(new JLabel("Sky Outer Radius:"));
        JSpinner jSpinner3 = new JSpinner(new SpinnerNumberModel(radialProfile.skyOuterRadius, 1, 100, 1));
        jPanel.add(jSpinner3);
        JButton jButton = new JButton("Update Plot");
        jButton.addActionListener(actionEvent -> {
            radialProfile.objectRadius = (Integer)jSpinner.getValue();
            radialProfile.skyInnerRadius = (Integer)jSpinner2.getValue();
            radialProfile.skyOuterRadius = (Integer)jSpinner3.getValue();
            radialProfilePlot.repaint();
        });
        jPanel.add(new JLabel());
        jPanel.add(jButton);
        jDialog.add((Component)jPanel, "South");
        jDialog.pack();
        jDialog.setLocationRelativeTo(this);
        jDialog.setVisible(true);
    }

    private TargetInfo fetchVSXTargetInfo(String string) {
        try {
            String string2;
            String string3 = String.format("https://vsx.aavso.org/index.php?view=api.object&ident=%s", URLEncoder.encode(string, "UTF-8"));
            URL uRL = new URL(string3);
            HttpURLConnection httpURLConnection = (HttpURLConnection)uRL.openConnection();
            httpURLConnection.setRequestMethod("GET");
            httpURLConnection.setConnectTimeout(10000);
            httpURLConnection.setReadTimeout(10000);
            int n = httpURLConnection.getResponseCode();
            if (n != 200) {
                System.err.println("VSX returned HTTP " + n);
                return null;
            }
            BufferedReader bufferedReader = new BufferedReader(new InputStreamReader(httpURLConnection.getInputStream()));
            StringBuilder stringBuilder = new StringBuilder();
            while ((string2 = bufferedReader.readLine()) != null) {
                stringBuilder.append(string2).append("\n");
            }
            bufferedReader.close();
            return this.parseVSXResponse(stringBuilder.toString());
        }
        catch (Exception exception) {
            System.err.println("Error fetching VSX data: " + exception.getMessage());
            return null;
        }
    }

    private TargetInfo parseVSXResponse(String string) {
        try {
            String string2;
            TargetInfo targetInfo = new TargetInfo("");
            String string3 = this.extractXMLValue(string, "Name");
            if (string3 == null || string3.isEmpty()) {
                return null;
            }
            targetInfo.name = string3;
            targetInfo.auid = this.extractXMLValue(string, "AUID");
            String string4 = this.extractXMLValue(string, "RA2000");
            String string5 = this.extractXMLValue(string, "Declination2000");
            if (string4 != null) {
                try {
                    targetInfo.ra = Double.parseDouble(string4);
                }
                catch (NumberFormatException numberFormatException) {
                    // empty catch block
                }
            }
            if (string5 != null) {
                try {
                    targetInfo.dec = Double.parseDouble(string5);
                }
                catch (NumberFormatException numberFormatException) {
                    // empty catch block
                }
            }
            targetInfo.variabilityType = this.extractXMLValue(string, "VariabilityType");
            String string6 = this.extractXMLValue(string, "Period");
            if (string6 != null) {
                try {
                    targetInfo.period = Double.parseDouble(string6);
                }
                catch (NumberFormatException numberFormatException) {
                    // empty catch block
                }
            }
            if ((string2 = this.extractXMLValue(string, "Epoch")) != null) {
                try {
                    targetInfo.epoch = Double.parseDouble(string2);
                }
                catch (NumberFormatException numberFormatException) {
                    // empty catch block
                }
            }
            targetInfo.maxMag = this.extractXMLValue(string, "MaxMag");
            targetInfo.minMag = this.extractXMLValue(string, "MinMag");
            targetInfo.constellation = this.extractXMLValue(string, "Constellation");
            return targetInfo;
        }
        catch (Exception exception) {
            System.err.println("Error parsing VSX XML: " + exception.getMessage());
            return null;
        }
    }

    private String extractXMLValue(String string, String string2) {
        String string3 = "<" + string2 + ">";
        String string4 = "</" + string2 + ">";
        int n = string.indexOf(string3);
        if (n < 0) {
            return null;
        }
        int n2 = string.indexOf(string4, n += string3.length());
        if (n2 < 0) {
            return null;
        }
        return string.substring(n, n2).trim();
    }

    private void fetchVSXFieldStars(double d, double d2, double d3) {
        this.variableStarsInField.clear();
        try {
            String string;
            int n = (int)(d / 15.0);
            double d4 = (d / 15.0 - (double)n) * 60.0;
            int n2 = (int)d4;
            double d5 = (d4 - (double)n2) * 60.0;
            String string2 = d2 >= 0.0 ? "+" : "-";
            double d6 = Math.abs(d2);
            int n3 = (int)d6;
            double d7 = (d6 - (double)n3) * 60.0;
            int n4 = (int)d7;
            double d8 = (d7 - (double)n4) * 60.0;
            String string3 = String.format("%02d %02d %05.2f %s%02d %02d %04.1f", n, n2, d5, string2, n3, n4, d8);
            String string4 = "ql=2&getCoordinates=0&plotType=Search&special=" + URLEncoder.encode("index.php?view=results.special&sid=2", "UTF-8") + "&ident=&constid=0&targetcenter=" + URLEncoder.encode(string3, "UTF-8") + "&format=s&fieldsize=" + (int)d3 + "&fieldunit=2&geometry=b&maxlower=&maxupper=&minlower=&minupper=&periodlower=&periodupper=&epochlower=&epochupper=&riselower=&riseupper=&novalower=&novaupper=&vartype=&vargroup=0&spectype=&assoctype=0&campaign=0&color=0&colorlower=&colorupper=&filter%5B%5D=0&filter%5B%5D=1&filter%5B%5D=2&filter%5B%5D=3&order=9";
            String string5 = "https://vsx.aavso.org/index.php?view=results.submit2";
            System.out.println("=== VSX Field Star Query ===");
            System.out.println("URL: " + string5);
            System.out.println("Center coordinates: " + string3 + " (RA=" + d + "\u00b0, Dec=" + d2 + "\u00b0)");
            System.out.println("FOV: " + d3 + " arcmin");
            System.out.println("Mag Limit: " + this.currentMagLimit);
            System.out.println("POST data: " + string4);
            URL uRL = new URL(string5);
            HttpURLConnection httpURLConnection = (HttpURLConnection)uRL.openConnection();
            httpURLConnection.setRequestMethod("POST");
            httpURLConnection.setDoOutput(true);
            httpURLConnection.setRequestProperty("Content-Type", "application/x-www-form-urlencoded");
            httpURLConnection.setRequestProperty("Content-Length", String.valueOf(string4.length()));
            httpURLConnection.setRequestProperty("User-Agent", "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7)");
            OutputStream outputStream = httpURLConnection.getOutputStream();
            outputStream.write(string4.getBytes("UTF-8"));
            outputStream.flush();
            outputStream.close();
            int n5 = httpURLConnection.getResponseCode();
            System.out.println("Response Code: " + n5);
            System.out.println("Response Message: " + httpURLConnection.getResponseMessage());
            BufferedReader bufferedReader = new BufferedReader(new InputStreamReader(httpURLConnection.getInputStream()));
            StringBuilder stringBuilder = new StringBuilder();
            int n6 = 0;
            while ((string = bufferedReader.readLine()) != null) {
                stringBuilder.append(string).append("\n");
                ++n6;
            }
            bufferedReader.close();
            System.out.println("Response length: " + stringBuilder.length() + " characters, " + n6 + " lines");
            try {
                FileWriter fileWriter = new FileWriter("/tmp/vsx_response.html");
                fileWriter.write(stringBuilder.toString());
                fileWriter.close();
                System.out.println("Response saved to: /tmp/vsx_response.html");
            }
            catch (Exception exception) {
                System.err.println("Could not save response: " + exception.getMessage());
            }
            this.parseVSXFieldResults(stringBuilder.toString());
            this.filterVariablesOnImage();
            System.out.println("Parsing complete. Found " + this.variableStarsInField.size() + " variable stars in field");
            System.out.println("=== End VSX Query ===");
        }
        catch (Exception exception) {
            System.err.println("Error fetching VSX field stars: " + exception.getMessage());
            exception.printStackTrace();
        }
    }

    private void filterVariablesOnImage() {
        if (this.variableStarsInField == null || this.variableStarsInField.isEmpty() || this.currentFitsHeader == null) {
            return;
        }
        try {
            Integer n = this.currentFitsHeader.getIntValue("NAXIS1");
            Integer n2 = this.currentFitsHeader.getIntValue("NAXIS2");
            if (n == null || n2 == null) {
                return;
            }
            ArrayList<TargetInfo> arrayList = new ArrayList<TargetInfo>();
            double d = Double.MAX_VALUE;
            for (TargetInfo targetInfo : this.variableStarsInField) {
                int[] nArray = this.raDecToPixel(targetInfo.ra, targetInfo.dec);
                if (nArray == null || nArray[0] < 0 || nArray[0] >= n || nArray[1] < 0 || nArray[1] >= n2) continue;
                arrayList.add(targetInfo);
                if (targetInfo.maxMag == null) continue;
                try {
                    double d2 = Double.parseDouble(targetInfo.maxMag.replaceAll("[^0-9.-]", ""));
                    if (!(d2 < d)) continue;
                    d = d2;
                }
                catch (Exception exception) {}
            }
            int n3 = this.variableStarsInField.size() - arrayList.size();
            this.variableStarsInField = arrayList;
            if (n3 > 0) {
                System.out.println("Filtered out " + n3 + " variables not on image");
            }
            if (d < Double.MAX_VALUE) {
                this.currentMagLimit = d + 2.0;
                System.out.println("Auto mag limit: " + String.format("%.2f", this.currentMagLimit) + " (brightest var max mag " + String.format("%.2f", d) + " + 2)");
                SwingUtilities.invokeLater(() -> this.magLimitSpinner.setValue(this.currentMagLimit));
            }
        }
        catch (Exception exception) {
            System.err.println("Error filtering variables on image: " + exception.getMessage());
        }
    }

    private void parseVSXFieldResults(String string) {
        System.out.println("=== Parsing VSX Results ===");
        String[] stringArray = string.split("\n");
        System.out.println("Total lines in HTML: " + stringArray.length);
        boolean bl = false;
        boolean bl2 = false;
        String string2 = null;
        String string3 = null;
        String string4 = null;
        String string5 = null;
        int n = 0;
        int n2 = 0;
        int n3 = 0;
        StringBuilder stringBuilder = new StringBuilder();
        for (int i = 0; i < stringArray.length; ++i) {
            String string6 = stringArray[i].trim();
            if (string6.startsWith("<tr>") || string6.startsWith("<tr style=")) {
                ++n3;
                bl = true;
                string2 = null;
                string3 = null;
                string4 = null;
                string5 = null;
                n = 0;
                if (++n2 < 19 || n2 > 23) continue;
                System.out.println("Found <tr> at line " + i + ", rowCount=" + n2);
                continue;
            }
            if (bl && n2 >= 19 && n2 <= 23 && i < stringArray.length && string6.startsWith("<td")) {
                System.out.println("  Line " + i + " in row " + n2 + ": [" + string6.substring(0, Math.min(60, string6.length())) + "...]");
            }
            if (string6.equals("</tr>") && bl) {
                if (n2 >= 19 && n2 <= 23) {
                    System.out.println("End of row " + n2 + ": name=" + string2 + ", coords=" + string3 + ", type=" + string4 + ", mag=" + string5 + ", cellIndex=" + n);
                }
                if (string2 != null && string3 != null && string5 != null) {
                    System.out.println("Row " + n2 + ": name=" + string2 + ", coords=" + string3 + ", type=" + string4 + ", mag=" + string5);
                    double[] dArray = this.parseCoordinates(string3);
                    if (dArray != null) {
                        System.out.println("  Parsed coords: RA=" + dArray[0] + "\u00b0, Dec=" + dArray[1] + "\u00b0");
                        double d = this.parseBrightestMagnitude(string5);
                        System.out.println("  Brightest mag: " + d + ", limit: " + this.currentMagLimit);
                        if (d <= this.currentMagLimit) {
                            String[] stringArray2;
                            TargetInfo targetInfo = new TargetInfo(string2);
                            targetInfo.ra = dArray[0];
                            targetInfo.dec = dArray[1];
                            String string7 = targetInfo.variabilityType = string4 != null ? string4 : "";
                            if (string5 != null && string5.contains("-") && (stringArray2 = string5.split("-")).length >= 2) {
                                String string8;
                                targetInfo.minMag = stringArray2[0].trim();
                                targetInfo.maxMag = string8 = stringArray2[1].trim().replaceAll("[A-Za-z]", "").trim();
                            }
                            this.variableStarsInField.add(targetInfo);
                            System.out.println("  \u2713 Added variable: " + string2 + " (mag " + d + ")");
                        } else {
                            System.out.println("  \u2717 Skipped (too faint): " + string2 + " (mag " + d + " > limit " + this.currentMagLimit + ")");
                        }
                    } else {
                        System.out.println("  \u2717 Failed to parse coordinates: " + string3);
                    }
                }
                bl = false;
                continue;
            }
            if (bl && string6.startsWith("<td class=\"indexdata\"")) {
                ++n;
                stringBuilder = new StringBuilder();
                if (n2 >= 19 && n2 <= 23) {
                    System.out.println("  Cell " + n + " started at line " + i);
                }
                if (string6.contains("</td>")) {
                    int n4 = string6.indexOf(">") + 1;
                    int n5 = string6.indexOf("</td>");
                    if (n5 > n4) {
                        String string9 = string6.substring(n4, n5).trim();
                        if (n2 >= 19 && n2 <= 23) {
                            System.out.println("  Cell " + n + " content (single-line): [" + string9 + "]");
                        }
                        if (n == 3 && string9.contains("<a href=")) {
                            int n6;
                            int n7 = string9.indexOf("<a href=");
                            if (n7 >= 0 && (n6 = string9.indexOf("</a>", n7)) > n7) {
                                int n8 = string9.indexOf(">", n7) + 1;
                                string2 = string9.substring(n8, n6).trim();
                                if (n2 >= 19 && n2 <= 23) {
                                    System.out.println("    \u2192 Extracted name: " + string2);
                                }
                            }
                        } else if (n == 5) {
                            string3 = string9.replaceAll("<[^>]*>", "").trim();
                            if (n2 >= 19 && n2 <= 23) {
                                System.out.println("    \u2192 Extracted coords: " + string3);
                            }
                        } else if (n == 7) {
                            string4 = string9.replaceAll("<[^>]*>", "").trim();
                            if (n2 >= 19 && n2 <= 23) {
                                System.out.println("    \u2192 Extracted type: " + string4);
                            }
                        } else if (n == 9) {
                            int n9 = string9.indexOf("<br>");
                            if (n9 >= 0) {
                                string9 = string9.substring(0, n9);
                            }
                            string5 = string9.replaceAll("<[^>]*>", "").trim();
                            if (n2 >= 19 && n2 <= 23) {
                                System.out.println("    \u2192 Extracted mag: " + string5);
                            }
                        }
                    }
                    bl2 = false;
                    continue;
                }
                bl2 = true;
                int n10 = string6.indexOf(">") + 1;
                if (n10 >= string6.length() || string6.substring(n10).trim().isEmpty()) continue;
                stringBuilder.append(string6.substring(n10));
                continue;
            }
            if (bl2 && string6.contains("</td>")) {
                String string10;
                bl2 = false;
                int n11 = string6.indexOf("</td>");
                if (n11 > 0 && !(string10 = string6.substring(0, n11).trim()).isEmpty()) {
                    if (stringBuilder.length() > 0) {
                        stringBuilder.append(" ");
                    }
                    stringBuilder.append(string10);
                }
                String string11 = stringBuilder.toString().trim();
                if (n2 >= 19 && n2 <= 23 && n <= 9) {
                    System.out.println("  Cell " + n + " content: [" + string11 + "]");
                }
                if (n == 3 && string11.contains("<a href=")) {
                    int n12;
                    int n13 = string11.indexOf("<a href=");
                    if (n13 < 0 || (n12 = string11.indexOf("</a>", n13)) <= n13) continue;
                    int n14 = string11.indexOf(">", n13) + 1;
                    string2 = string11.substring(n14, n12).trim();
                    if (n2 < 19 || n2 > 23) continue;
                    System.out.println("    \u2192 Extracted name: " + string2);
                    continue;
                }
                if (n == 5) {
                    string3 = string11.replaceAll("<[^>]*>", "").trim();
                    if (n2 < 19 || n2 > 23) continue;
                    System.out.println("    \u2192 Extracted coords: " + string3);
                    continue;
                }
                if (n == 7) {
                    string4 = string11.replaceAll("<[^>]*>", "").trim();
                    if (n2 < 19 || n2 > 23) continue;
                    System.out.println("    \u2192 Extracted type: " + string4);
                    continue;
                }
                if (n != 9) continue;
                int n15 = string11.indexOf("<br>");
                if (n15 >= 0) {
                    string11 = string11.substring(0, n15);
                }
                string5 = string11.replaceAll("<[^>]*>", "").trim();
                if (n2 < 19 || n2 > 23) continue;
                System.out.println("    \u2192 Extracted mag: " + string5);
                continue;
            }
            if (!bl2 || string6.isEmpty()) continue;
            if (stringBuilder.length() > 0) {
                stringBuilder.append(" ");
            }
            stringBuilder.append(string6);
        }
        System.out.println("Total <tr> rows found: " + n3);
        System.out.println("Total data rows processed: " + n2);
        System.out.println("Variables added: " + this.variableStarsInField.size());
    }

    private double[] parseCoordinates(String string) {
        try {
            String[] stringArray = string.trim().split("\\s+");
            if (stringArray.length >= 6) {
                double d = Double.parseDouble(stringArray[0]);
                double d2 = Double.parseDouble(stringArray[1]);
                double d3 = Double.parseDouble(stringArray[2]);
                double d4 = (d + d2 / 60.0 + d3 / 3600.0) * 15.0;
                String string2 = stringArray[3];
                boolean bl = string2.startsWith("-");
                double d5 = Math.abs(Double.parseDouble(string2));
                double d6 = Double.parseDouble(stringArray[4]);
                double d7 = Double.parseDouble(stringArray[5]);
                double d8 = d5 + d6 / 60.0 + d7 / 3600.0;
                if (bl) {
                    d8 = -d8;
                }
                return new double[]{d4, d8};
            }
        }
        catch (Exception exception) {
            System.err.println("Failed to parse coordinates: " + string);
        }
        return null;
    }

    private double parseBrightestMagnitude(String string) {
        try {
            String[] stringArray = string.split("-");
            if (stringArray.length >= 1) {
                String string2 = stringArray[0].trim();
                if (!(string2 = string2.replaceAll("[^0-9.]", "")).isEmpty()) {
                    return Double.parseDouble(string2);
                }
            }
        }
        catch (Exception exception) {
            System.err.println("Failed to parse magnitude range: " + string);
        }
        return 99.9;
    }

    private List<ReferenceStar> fetchAAVSOVSD(String string, Double d, Double d2, double d3, double d4) {
        try {
            String string2;
            Object object;
            Object object2;
            String string3;
            if (string != null && !string.isEmpty()) {
                string3 = String.format("https://apps.aavso.org/vsp/photometry/?star=%s&fov=%.1f&maglimit=%.1f", URLEncoder.encode(string, "UTF-8"), d3, d4);
            } else if (d != null && d2 != null) {
                object2 = this.formatRA(d);
                object = this.formatDec(d2);
                string3 = String.format("https://apps.aavso.org/vsp/photometry/?ra=%s&dec=%s&fov=%.1f&maglimit=%.1f", URLEncoder.encode((String)object2, "UTF-8"), URLEncoder.encode((String)object, "UTF-8"), d3, d4);
                System.out.println("Using coordinates - URL: " + string3);
            } else {
                System.err.println("Either targetName or ra/dec must be provided");
                return null;
            }
            object2 = new URL(string3);
            object = (HttpURLConnection)((URL)object2).openConnection();
            ((HttpURLConnection)object).setRequestMethod("GET");
            ((URLConnection)object).setConnectTimeout(10000);
            ((URLConnection)object).setReadTimeout(10000);
            int n = ((HttpURLConnection)object).getResponseCode();
            if (n != 200) {
                System.err.println("AAVSO VSD returned HTTP " + n);
                return null;
            }
            BufferedReader bufferedReader = new BufferedReader(new InputStreamReader(((URLConnection)object).getInputStream()));
            StringBuilder stringBuilder = new StringBuilder();
            while ((string2 = bufferedReader.readLine()) != null) {
                stringBuilder.append(string2).append("\n");
            }
            bufferedReader.close();
            return this.parseAAVSOVSD(stringBuilder.toString());
        }
        catch (Exception exception) {
            System.err.println("Error fetching AAVSO VSD: " + exception.getMessage());
            exception.printStackTrace();
            return null;
        }
    }

    private List<ReferenceStar> parseAAVSOVSD(String string) {
        ArrayList<ReferenceStar> arrayList = new ArrayList<ReferenceStar>();
        try {
            String[] stringArray = string.split("\n");
            boolean bl = false;
            boolean bl2 = false;
            String string2 = null;
            String string3 = null;
            String string4 = null;
            String string5 = null;
            String string6 = null;
            String string7 = null;
            int n = 0;
            for (String string8 : stringArray) {
                if ((string8 = string8.trim()).contains("Source Reference Table")) break;
                if (string8.contains("<table")) {
                    bl = true;
                    continue;
                }
                if (string8.contains("</table>")) {
                    bl = false;
                    continue;
                }
                if (!bl) continue;
                if (string8.contains("<tr>")) {
                    bl2 = true;
                    n = 0;
                    string2 = null;
                    string3 = null;
                    string4 = null;
                    string5 = null;
                    string6 = null;
                    string7 = null;
                    continue;
                }
                if (string8.contains("</tr>")) {
                    bl2 = false;
                    if (string2 == null || string3 == null || string4 == null || !string2.matches(".*\\d.*") || !string2.contains("-")) continue;
                    try {
                        double d = this.extractDecimalDegrees(string3);
                        double d2 = this.extractDecimalDegrees(string4);
                        double d3 = 0.0;
                        double d4 = 0.0;
                        if (string6 != null) {
                            double[] dArray = this.parseMagWithError(string6);
                            d3 = dArray[0];
                            d4 = dArray[1];
                        }
                        double d5 = 0.0;
                        double d6 = 0.0;
                        if (string7 != null) {
                            double[] dArray = this.parseMagWithError(string7);
                            d5 = dArray[0];
                            d6 = dArray[1];
                        }
                        arrayList.add(new ReferenceStar(string2, d, d2, d3, d4, d5, d6, string5));
                    }
                    catch (Exception exception) {}
                    continue;
                }
                if (!bl2 || !string8.contains("<td>")) continue;
                String string9 = string8.replaceAll("<[^>]*>", "").trim();
                switch (n) {
                    case 0: {
                        if (string9.isEmpty() || string9.equals("AUID")) break;
                        string2 = string9;
                        break;
                    }
                    case 1: {
                        if (string9.isEmpty() || string9.equals("RA")) break;
                        string3 = string9;
                        break;
                    }
                    case 2: {
                        if (string9.isEmpty() || string9.equals("Dec")) break;
                        string4 = string9;
                        break;
                    }
                    case 3: {
                        if (string9.isEmpty() || string9.equals("Label")) break;
                        string5 = string9;
                        break;
                    }
                    case 4: {
                        if (string9.isEmpty() || string9.equals("V")) break;
                        string6 = string9;
                        break;
                    }
                    case 5: {
                        if (string9.isEmpty() || string9.equals("B-V")) break;
                        string7 = string9;
                    }
                }
                ++n;
            }
            System.out.println("Parsed " + arrayList.size() + " comparison stars from AAVSO VSD");
            for (ReferenceStar referenceStar : arrayList) {
                System.out.println("  " + referenceStar.name + " (label=" + referenceStar.label + "): V=" + String.format("%.3f\u00b1%.3f", referenceStar.mag, referenceStar.magError) + ", B-V=" + String.format("%.3f\u00b1%.3f", referenceStar.bv, referenceStar.bvError));
            }
        }
        catch (Exception exception) {
            System.err.println("Error parsing AAVSO VSD HTML: " + exception.getMessage());
        }
        return arrayList;
    }

    private double[] parseMagWithError(String string) {
        double d = 0.0;
        double d2 = 0.0;
        try {
            string = string.replaceAll("<sup>.*?</sup>", "").trim();
            int n = string.indexOf(40);
            int n2 = string.indexOf(41);
            if (n >= 0 && n2 > n) {
                String string2 = string.substring(0, n).trim();
                d = Double.parseDouble(string2);
                String string3 = string.substring(n + 1, n2).trim();
                d2 = Double.parseDouble(string3);
            } else {
                d = Double.parseDouble(string.trim());
            }
        }
        catch (Exception exception) {
            // empty catch block
        }
        return new double[]{d, d2};
    }

    private double extractDecimalDegrees(String string) {
        int n = string.indexOf(91);
        int n2 = string.indexOf(93);
        if (n >= 0 && n2 > n) {
            String string2 = string.substring(n + 1, n2);
            string2 = string2.replace("\u00b0", "").replace("&deg;", "").trim();
            return Double.parseDouble(string2);
        }
        return this.parseRADec(string, string.contains(":") && !string.contains("+") && !string.contains("-"));
    }

    private double parseRADec(String string, boolean bl) {
        try {
            return Double.parseDouble(string);
        }
        catch (NumberFormatException numberFormatException) {
            String[] stringArray;
            String[] stringArray2 = stringArray = string.contains(":") ? string.split(":") : string.split("\\s+");
            if (stringArray.length == 3) {
                double d = Double.parseDouble(stringArray[0].trim());
                double d2 = Double.parseDouble(stringArray[1].trim());
                double d3 = Double.parseDouble(stringArray[2].trim());
                if (bl) {
                    return (d + d2 / 60.0 + d3 / 3600.0) * 15.0;
                }
                double d4 = d >= 0.0 ? 1.0 : -1.0;
                return d4 * (Math.abs(d) + d2 / 60.0 + d3 / 3600.0);
            }
            return 0.0;
        }
    }

    private Double getCenterRA(Header header) {
        try {
            Double d = header.getDoubleValue("CRVAL1");
            if (d != null && Math.abs(d) > 0.001) {
                System.out.println("Found CRVAL1 (RA) = " + d);
                return d;
            }
            d = header.getDoubleValue("RA");
            if (d != null && Math.abs(d) > 0.001) {
                System.out.println("Found RA = " + d);
                return d;
            }
            String string = header.getStringValue("OBJCTRA");
            if (string != null && !string.trim().isEmpty()) {
                d = this.parseRADec(string.trim(), true);
                System.out.println("Parsed OBJCTRA string '" + string + "' = " + d);
                return d;
            }
            string = header.getStringValue("RA");
            if (string != null && !string.trim().isEmpty()) {
                d = this.parseRADec(string.trim(), true);
                System.out.println("Parsed RA string = " + d);
                return d;
            }
            d = header.getDoubleValue("OBJCTRA");
            if (d != null && Math.abs(d) > 0.001) {
                System.out.println("Found OBJCTRA (RA) = " + d);
                return d;
            }
            System.err.println("Could not find RA in FITS header. Available keys:");
            header.iterator().forEachRemaining(headerCard -> {
                String string = headerCard.getKey();
                if (string.contains("RA") || string.contains("CRVAL") || string.equals("OBJCTRA")) {
                    System.err.println("  " + string + " = " + headerCard.getValue());
                }
            });
        }
        catch (Exception exception) {
            System.err.println("Error reading RA from header: " + exception.getMessage());
        }
        return null;
    }

    private Double getCenterDec(Header header) {
        try {
            Double d = header.getDoubleValue("CRVAL2");
            if (d != null && Math.abs(d) > 0.001) {
                System.out.println("Found CRVAL2 (Dec) = " + d);
                return d;
            }
            d = header.getDoubleValue("DEC");
            if (d != null && Math.abs(d) > 0.001) {
                System.out.println("Found DEC = " + d);
                return d;
            }
            String string = header.getStringValue("OBJCTDEC");
            if (string != null && !string.trim().isEmpty()) {
                d = this.parseRADec(string.trim(), false);
                System.out.println("Parsed OBJCTDEC string '" + string + "' = " + d);
                return d;
            }
            string = header.getStringValue("DEC");
            if (string != null && !string.trim().isEmpty()) {
                d = this.parseRADec(string.trim(), false);
                System.out.println("Parsed DEC string = " + d);
                return d;
            }
            d = header.getDoubleValue("OBJCTDEC");
            if (d != null && Math.abs(d) > 0.001) {
                System.out.println("Found OBJCTDEC (Dec) = " + d);
                return d;
            }
            System.err.println("Could not find DEC in FITS header. Available keys:");
            header.iterator().forEachRemaining(headerCard -> {
                String string = headerCard.getKey();
                if (string.contains("DEC") || string.contains("CRVAL") || string.equals("OBJCTDEC")) {
                    System.err.println("  " + string + " = " + headerCard.getValue());
                }
            });
        }
        catch (Exception exception) {
            System.err.println("Error reading Dec from header: " + exception.getMessage());
        }
        return null;
    }

    private String formatRA(double d) {
        double d2 = d / 15.0;
        int n = (int)d2;
        double d3 = (d2 - (double)n) * 60.0;
        int n2 = (int)d3;
        double d4 = (d3 - (double)n2) * 60.0;
        return String.format("%02d:%02d:%05.2f", n, n2, d4);
    }

    private String formatDec(double d) {
        String string = d >= 0.0 ? "+" : "-";
        double d2 = Math.abs(d);
        int n = (int)d2;
        double d3 = (d2 - (double)n) * 60.0;
        int n2 = (int)d3;
        double d4 = (d3 - (double)n2) * 60.0;
        return String.format("%s%02d:%02d:%05.2f", string, n, n2, d4);
    }

    private void resetZoom() {
        this.zoomLevel = 1.0;
        this.panOffsetX = 0;
        this.panOffsetY = 0;
        this.imagePanel.repaint();
    }

    private void performPhotometry() {
        if (this.currentImageData == null || this.currentFitsHeader == null) {
            JOptionPane.showMessageDialog(this, "No image loaded", "Cannot Perform Photometry", 2);
            return;
        }
        if (this.variableStarsInField == null || this.variableStarsInField.isEmpty()) {
            JOptionPane.showMessageDialog(this, "No variable stars loaded. Please wait for VSD query to complete.", "Cannot Perform Photometry", 2);
            return;
        }
        if (this.currentRadialProfile == null) {
            JOptionPane.showMessageDialog(this, "No radial profile available. Please load an image first.", "Cannot Perform Photometry", 2);
            return;
        }
        this.statusLabel.setText("Performing photometry...");
        this.statusLabel.setForeground(Color.BLUE);
        SwingWorker<List<PhotometryResult>, String> swingWorker = new SwingWorker<List<PhotometryResult>, String>(){

            @Override
            protected List<PhotometryResult> doInBackground() throws Exception {
                Object object;
                double d;
                ArrayList<PhotometryResult> arrayList = new ArrayList<PhotometryResult>();
                String string = TargetPhotometryPanel.this.currentFitsHeader.getStringValue("FILTER");
                if (string == null || string.trim().isEmpty()) {
                    string = "Unknown";
                }
                int n = TargetPhotometryPanel.this.currentRadialProfile.objectRadius;
                int n2 = TargetPhotometryPanel.this.currentRadialProfile.skyInnerRadius;
                int n3 = TargetPhotometryPanel.this.currentRadialProfile.skyOuterRadius;
                this.publish("Using aperture: obj=" + n + ", sky=" + n2 + "-" + n3);
                int n4 = 0;
                for (TargetInfo object2 : TargetPhotometryPanel.this.variableStarsInField) {
                    int[] referenceStar;
                    if (object2.maxMag != null && object2.minMag != null) {
                        try {
                            double exception = Double.parseDouble(object2.maxMag.replaceAll("[^0-9.-]", ""));
                            double n6 = Double.parseDouble(object2.minMag.replaceAll("[^0-9.-]", ""));
                            double fluxMeasurement = Math.abs(exception - n6);
                            if (fluxMeasurement < TargetPhotometryPanel.this.minVarDeltaV) {
                                this.publish("Skipping " + object2.name + " - amplitude " + String.format("%.3f", fluxMeasurement) + " < " + String.format("%.3f", TargetPhotometryPanel.this.minVarDeltaV));
                                continue;
                            }
                        }
                        catch (Exception nArray) {
                            // empty catch block
                        }
                    }
                    if ((referenceStar = TargetPhotometryPanel.this.raDecToPixel(object2.ra, object2.dec)) == null) {
                        this.publish("Skipping " + object2.name + " - could not convert coordinates");
                        continue;
                    }
                    int string3 = referenceStar[0];
                    int nArray = referenceStar[1];
                    if (string3 < n3 || nArray < n3 || string3 >= TargetPhotometryPanel.this.currentImageData.length - n3 || nArray >= TargetPhotometryPanel.this.currentImageData[0].length - n3) {
                        this.publish("Skipping " + object2.name + " - too close to edge");
                        continue;
                    }
                    double[] n8 = TargetPhotometryPanel.this.refineCentroid(TargetPhotometryPanel.this.currentImageData, string3, nArray, 10);
                    if (n8 == null) {
                        this.publish("Skipping " + object2.name + " - centroid refinement failed");
                        continue;
                    }
                    FluxMeasurement n9 = TargetPhotometryPanel.this.measureFlux(TargetPhotometryPanel.this.currentImageData, n8[0], n8[1], n, n2, n3, object2.name);
                    if (n9 == null) {
                        this.publish("Skipping " + object2.name + " - flux measurement failed");
                        continue;
                    }
                    double dArray = -2.5 * Math.log10(n9.netFlux);
                    d = 1.0857 * n9.error / n9.netFlux;
                    double[] d6 = TargetPhotometryPanel.this.pixelToRaDec(n8[0], n8[1]);
                    String string2 = TargetPhotometryPanel.this.currentFitsFile != null ? TargetPhotometryPanel.this.currentFitsFile.getName() : "unknown";
                    object = new PhotometryResult(string2, object2.name, "Variable", string, object2.ra, object2.dec, n8[0], n8[1], d6 != null ? d6[0] : object2.ra, d6 != null ? d6[1] : object2.dec, n9.starFlux, n9.skyFlux, n9.skyPerPixel, n9.netFlux, dArray, d);
                    arrayList.add((PhotometryResult)object);
                    this.publish("Processed variable star " + ++n4 + ": " + object2.name);
                }
                TargetPhotometryPanel.this.selectCheckStar();
                if (TargetPhotometryPanel.this.referenceStars != null) {
                    this.publish("Processing " + TargetPhotometryPanel.this.referenceStars.size() + " reference stars...");
                    int n7 = 0;
                    for (ReferenceStar referenceStar : TargetPhotometryPanel.this.referenceStars) {
                        Object object2;
                        FluxMeasurement n10;
                        String string3 = referenceStar.label.isEmpty() ? referenceStar.name : referenceStar.label;
                        int[] nArray = TargetPhotometryPanel.this.raDecToPixel(referenceStar.ra, referenceStar.dec);
                        if (nArray == null) {
                            this.publish("Skipping comp " + string3 + " - could not convert coordinates");
                            continue;
                        }
                        int n5 = nArray[0];
                        int n6 = nArray[1];
                        if (string3.equals("123") || string3.equals("129")) {
                            System.out.println("DEBUG " + string3 + ": initial pixel coords: x=" + n5 + ", y=" + n6);
                        }
                        if (n5 < n3 || n6 < n3 || n5 >= TargetPhotometryPanel.this.currentImageData.length - n3 || n6 >= TargetPhotometryPanel.this.currentImageData[0].length - n3) {
                            this.publish("Skipping comp " + string3 + " - too close to edge (x=" + n5 + ", y=" + n6 + ")");
                            continue;
                        }
                        double[] dArray = TargetPhotometryPanel.this.refineCentroid(TargetPhotometryPanel.this.currentImageData, n5, n6, 10);
                        if (dArray == null) {
                            this.publish("Skipping comp " + string3 + " - centroid refinement failed");
                            continue;
                        }
                        if (string3.equals("123") || string3.equals("129")) {
                            System.out.println("DEBUG " + string3 + ": refined centroid: x=" + dArray[0] + ", y=" + dArray[1]);
                            System.out.println("DEBUG " + string3 + ": shift from initial: dx=" + (dArray[0] - (double)n5) + ", dy=" + (dArray[1] - (double)n6));
                            int n8 = TargetPhotometryPanel.this.currentImageData[n5][n6];
                            System.out.println("DEBUG " + string3 + ": pixel value at initial position: " + n8);
                        }
                        if ((n10 = TargetPhotometryPanel.this.measureFlux(TargetPhotometryPanel.this.currentImageData, dArray[0], dArray[1], n, n2, n3, string3)) == null) {
                            this.publish("Skipping comp " + string3 + " - flux measurement failed");
                            continue;
                        }
                        d = -2.5 * Math.log10(n10.netFlux);
                        double d2 = 1.0857 * n10.error / n10.netFlux;
                        object = TargetPhotometryPanel.this.pixelToRaDec(dArray[0], dArray[1]);
                        String string4 = TargetPhotometryPanel.this.currentFitsFile != null ? TargetPhotometryPanel.this.currentFitsFile.getName() : "unknown";
                        String string5 = "Comparison";
                        if (TargetPhotometryPanel.this.checkStar != null) {
                            Object object3 = object2 = TargetPhotometryPanel.this.checkStar.label.isEmpty() ? TargetPhotometryPanel.this.checkStar.name : TargetPhotometryPanel.this.checkStar.label;
                            if (string3.equals(object2)) {
                                string5 = "Check";
                                System.out.println("CHECK STAR IDENTIFIED: " + string3 + " - setting type to 'Check'");
                            }
                        }
                        object2 = new PhotometryResult(string4, referenceStar.label.isEmpty() ? referenceStar.name : referenceStar.label, string5, string, referenceStar.ra, referenceStar.dec, dArray[0], dArray[1], object != null ? (double)object[0] : referenceStar.ra, object != null ? (double)object[1] : referenceStar.dec, n10.starFlux, n10.skyFlux, n10.skyPerPixel, n10.netFlux, d, d2);
                        arrayList.add((PhotometryResult)object2);
                        this.publish("Processed comp star " + ++n7 + ": " + string3);
                    }
                    this.publish("Total comp stars processed: " + n7 + " out of " + TargetPhotometryPanel.this.referenceStars.size());
                }
                this.publish("Photometry complete - " + arrayList.size() + " stars measured");
                if (TargetPhotometryPanel.this.referenceStars != null && !TargetPhotometryPanel.this.referenceStars.isEmpty()) {
                    TargetPhotometryPanel.this.calculateZeropoints(arrayList);
                }
                return arrayList;
            }

            @Override
            protected void process(List<String> list) {
                for (String string : list) {
                    System.out.println("Photometry: " + string);
                }
            }

            @Override
            protected void done() {
                try {
                    List list = (List)this.get();
                    if (TargetPhotometryPanel.this.photometryResults == null) {
                        TargetPhotometryPanel.this.photometryResults = new ArrayList<PhotometryResult>();
                    }
                    TargetPhotometryPanel.this.photometryResults.addAll(list);
                    TargetPhotometryPanel.this.updatePhotometryResultsTable();
                    TargetPhotometryPanel.this.statusLabel.setText("Photometry complete - " + TargetPhotometryPanel.this.photometryResults.size() + " total measurements (" + list.size() + " from this image)");
                    TargetPhotometryPanel.this.statusLabel.setForeground(new Color(0, 128, 0));
                }
                catch (Exception exception) {
                    exception.printStackTrace();
                    TargetPhotometryPanel.this.statusLabel.setText("Photometry failed: " + exception.getMessage());
                    TargetPhotometryPanel.this.statusLabel.setForeground(Color.RED);
                }
            }
        };
        swingWorker.execute();
    }

    private void selectCheckStar() {
        double d;
        double d2;
        Object object;
        if (this.referenceStars == null || this.referenceStars.isEmpty()) {
            System.out.println("No comparison stars available for CHECK star selection");
            this.checkStar = null;
            return;
        }
        if (this.currentImageData == null) {
            System.out.println("No image loaded for CHECK star selection");
            this.checkStar = null;
            return;
        }
        System.out.println("\n=== CHECK Star Selection ===");
        double d3 = 99.9;
        String string = "unknown";
        if (this.variableStarsInField != null && !this.variableStarsInField.isEmpty()) {
            object = null;
            d2 = 99.9;
            for (TargetInfo object2 : this.variableStarsInField) {
                if (object2.maxMag == null || object2.minMag == null) continue;
                try {
                    double d4;
                    double d5 = Double.parseDouble(object2.maxMag.replaceAll("[^0-9.-]", ""));
                    d = (d5 + (d4 = Double.parseDouble(object2.minMag.replaceAll("[^0-9.-]", "")))) / 2.0;
                    if (!(d < d2)) continue;
                    d2 = d;
                    object = object2;
                }
                catch (Exception exception) {}
            }
            if (object != null) {
                string = ((TargetInfo)object).name;
                d3 = d2;
            }
        }
        if (d3 >= 90.0) {
            System.out.println("No valid variable star magnitude available for CHECK star selection");
            System.out.println("Selecting brightest comparison star as CHECK");
            object = null;
            d2 = 99.9;
            for (ReferenceStar referenceStar : this.referenceStars) {
                if (!(referenceStar.mag < d2)) continue;
                d2 = referenceStar.mag;
                object = referenceStar;
            }
            this.checkStar = object;
            if (this.checkStar != null) {
                System.out.println("Selected CHECK star (brightest comp): " + (this.checkStar.label.isEmpty() ? this.checkStar.name : this.checkStar.label) + " (mag=" + String.format("%.3f", this.checkStar.mag) + ")");
            }
            System.out.println("=== End CHECK Star Selection ===\n");
            return;
        }
        System.out.println("Reference: " + string);
        if (d3 < 90.0) {
            System.out.println("Reference magnitude (avg): " + String.format("%.3f", d3));
        }
        int n = this.currentImageData.length;
        int n2 = this.currentImageData[0].length;
        double d6 = (double)n / 2.0;
        double d7 = (double)n2 / 2.0;
        double d8 = (double)Math.min(n, n2) / 2.0;
        System.out.println("Image dimensions: " + n + " x " + n2);
        System.out.println("Image center: (" + String.format("%.1f", d6) + ", " + String.format("%.1f", d7) + ")");
        System.out.println("Max distance from center: " + String.format("%.1f", d8) + " pixels");
        ReferenceStar referenceStar = null;
        d = Double.MAX_VALUE;
        for (ReferenceStar referenceStar2 : this.referenceStars) {
            int[] nArray = this.raDecToPixel(referenceStar2.ra, referenceStar2.dec);
            if (nArray == null) continue;
            double d9 = nArray[0];
            double d10 = nArray[1];
            double d11 = Math.sqrt(Math.pow(d9 - d6, 2.0) + Math.pow(d10 - d7, 2.0));
            if (d11 > d8) {
                System.out.println("  " + (referenceStar2.label.isEmpty() ? referenceStar2.name : referenceStar2.label) + ": mag=" + String.format("%.3f", referenceStar2.mag) + ", dist=" + String.format("%.1f", d11) + " > " + String.format("%.1f", d8) + " (too far)");
                continue;
            }
            double d12 = Math.abs(referenceStar2.mag - d3);
            System.out.println("  " + (referenceStar2.label.isEmpty() ? referenceStar2.name : referenceStar2.label) + ": mag=" + String.format("%.3f", referenceStar2.mag) + ", \u0394mag=" + String.format("%.3f", d12) + ", dist=" + String.format("%.1f", d11) + " pixels");
            if (!(d12 < d)) continue;
            d = d12;
            referenceStar = referenceStar2;
        }
        this.checkStar = referenceStar;
        if (this.checkStar != null) {
            System.out.println("\nSelected CHECK star: " + (this.checkStar.label.isEmpty() ? this.checkStar.name : this.checkStar.label) + " (mag=" + String.format("%.3f", this.checkStar.mag) + ", \u0394mag=" + String.format("%.3f", d) + ")");
        } else {
            System.out.println("\nNo suitable CHECK star found");
        }
        System.out.println("=== End CHECK Star Selection ===\n");
    }

    private void calculateZeropoints(List<PhotometryResult> list) {
        double d;
        if (this.referenceStars == null || this.referenceStars.isEmpty()) {
            System.out.println("No comparison stars for zeropoint calculation");
            return;
        }
        System.out.println("\n=== Zeropoint Calculation ===");
        String string = this.currentFitsFile != null ? this.currentFitsFile.getName() : "unknown";
        System.out.println("Image: " + string);
        ArrayList<Double> arrayList = new ArrayList<Double>();
        ArrayList<Double> arrayList2 = new ArrayList<Double>();
        ArrayList<String> arrayList3 = new ArrayList<String>();
        int n = 0;
        int n2 = 0;
        for (PhotometryResult photometryResult : list) {
            if (!photometryResult.type.equals("Comparison")) continue;
            ReferenceStar referenceStar = null;
            for (ReferenceStar referenceStar2 : this.referenceStars) {
                String string2 = referenceStar2.label.isEmpty() ? referenceStar2.name : referenceStar2.label;
                if (!string2.equals(photometryResult.name)) continue;
                referenceStar = referenceStar2;
                break;
            }
            if (referenceStar == null) {
                System.out.println("  WARNING: No reference star found for " + photometryResult.name);
                continue;
            }
            if (this.checkStar != null && referenceStar == this.checkStar) {
                System.out.println("  " + photometryResult.name + " (CHECK): excluded from zeropoint calculation");
                ++n2;
                continue;
            }
            double d2 = referenceStar.mag - photometryResult.instMag;
            double d3 = Math.max(photometryResult.magError, 0.001);
            d = 1.0 / (d3 * d3);
            arrayList.add(d2);
            arrayList2.add(d);
            arrayList3.add(photometryResult.name);
            ++n;
            System.out.println(String.format("  %s: std_mag=%.3f, inst_mag=%.3f, ZP=%.3f, err=%.4f, weight=%.1f", photometryResult.name, referenceStar.mag, photometryResult.instMag, d2, photometryResult.magError, d));
        }
        if (arrayList.isEmpty()) {
            System.out.println("No comparison stars available for zeropoint calculation");
            System.out.println("=== End Zeropoint Calculation ===\n");
            return;
        }
        double d4 = this.calculateWeightedMedian(arrayList, arrayList2);
        double d5 = 0.0;
        double d6 = 0.0;
        for (int i = 0; i < arrayList.size(); ++i) {
            d = (Double)arrayList.get(i) - d4;
            d6 += (Double)arrayList2.get(i) * d * d;
            d5 += ((Double)arrayList2.get(i)).doubleValue();
        }
        double d7 = Math.sqrt(d6 / d5);
        System.out.println("\n" + n + " comparison stars used (" + n2 + " excluded)");
        System.out.println(String.format("Weighted Median Zeropoint: %.4f \u00b1 %.4f", d4, d7));
        System.out.println("\n=== Applying Zeropoint Calibration ===");
        System.out.println("Variable Stars:");
        for (PhotometryResult photometryResult : list) {
            if (!photometryResult.type.equals("Variable")) continue;
            photometryResult.applyZeropoint(d4, d7);
            System.out.println(String.format("  %s: inst_mag=%.3f + ZP(%.4f) = calibrated_mag=%.3f \u00b1 %.4f", photometryResult.name, photometryResult.instMag, d4, photometryResult.calibratedMag, photometryResult.magError));
        }
        System.out.println("\nComparison Stars:");
        for (PhotometryResult photometryResult : list) {
            if (!photometryResult.type.equals("Comparison") && !photometryResult.type.equals("Check")) continue;
            ReferenceStar referenceStar = null;
            for (ReferenceStar referenceStar3 : this.referenceStars) {
                String string3 = referenceStar3.label.isEmpty() ? referenceStar3.name : referenceStar3.label;
                if (!string3.equals(photometryResult.name)) continue;
                referenceStar = referenceStar3;
                break;
            }
            if (referenceStar == null) continue;
            photometryResult.applyZeropoint(d4, d7);
            double d8 = photometryResult.calibratedMag - referenceStar.mag;
            if (photometryResult.type.equals("Check")) {
                System.out.println(String.format("  %s (CHECK): inst_mag=%.3f + ZP(%.4f) = calibrated_mag=%.3f \u00b1 %.4f", photometryResult.name, photometryResult.instMag, d4, photometryResult.calibratedMag, photometryResult.magError));
                System.out.println(String.format("    Standard mag=%.3f, Residual=%.4f (calibrated - standard)", referenceStar.mag, d8));
                if (!(Math.abs(d8) > 0.1)) continue;
                System.out.println("    WARNING: Large residual detected! CHECK star deviates by > 0.1 mag");
                continue;
            }
            System.out.println(String.format("  %s: inst_mag=%.3f + ZP(%.4f) = calibrated_mag=%.3f \u00b1 %.4f, residual=%.4f", photometryResult.name, photometryResult.instMag, d4, photometryResult.calibratedMag, photometryResult.magError, d8));
        }
        System.out.println("=== End Zeropoint Calibration ===\n");
    }

    private double calculateWeightedMedian(List<Double> list, List<Double> list2) {
        if (list.isEmpty()) {
            return 0.0;
        }
        if (list.size() == 1) {
            return list.get(0);
        }
        Integer[] integerArray = new Integer[list.size()];
        for (int i = 0; i < integerArray.length; ++i) {
            integerArray[i] = i;
        }
        Arrays.sort(integerArray, (n, n2) -> Double.compare((Double)list.get((int)n), (Double)list.get((int)n2)));
        double d = 0.0;
        for (double d2 : list2) {
            d += d2;
        }
        double d3 = d / 2.0;
        double d4 = 0.0;
        for (int i = 0; i < integerArray.length; ++i) {
            int n3 = integerArray[i];
            if (!((d4 += list2.get(n3).doubleValue()) >= d3)) continue;
            return list.get(n3);
        }
        return list.get(integerArray[integerArray.length / 2]);
    }

    private void performPhotometryAll() {
        File[] fileArray;
        if (this.currentTargetDir == null || !this.currentTargetDir.exists()) {
            JOptionPane.showMessageDialog(this, "No target directory loaded. Please select a target first.", "Cannot Perform Batch Photometry", 2);
            return;
        }
        if (this.variableStarsInField == null || this.variableStarsInField.isEmpty()) {
            JOptionPane.showMessageDialog(this, "No variable stars loaded. Please wait for VSD query to complete.", "Cannot Perform Batch Photometry", 2);
            return;
        }
        if (this.currentRadialProfile == null) {
            JOptionPane.showMessageDialog(this, "No radial profile available. Please load an image first.", "Cannot Perform Batch Photometry", 2);
            return;
        }
        final ArrayList<File> arrayList = new ArrayList<File>();
        File file2 = new File(this.currentTargetDir, "Reduced_images");
        if (file2.exists() && (fileArray = file2.listFiles((file, string) -> string.toLowerCase().endsWith("_cal.fits"))) != null && fileArray.length > 0) {
            arrayList.addAll(Arrays.asList(fileArray));
        }
        if (arrayList.isEmpty() && (fileArray = this.currentTargetDir.listFiles((file, string) -> string.toLowerCase().endsWith(".fits"))) != null && fileArray.length > 0) {
            arrayList.addAll(Arrays.asList(fileArray));
        }
        if (arrayList.isEmpty()) {
            JOptionPane.showMessageDialog(this, "No FITS files found in " + this.currentTargetDir.getName(), "Cannot Perform Batch Photometry", 2);
            return;
        }
        final int n = arrayList.size();
        System.out.println("========================================");
        System.out.println("BATCH PHOTOMETRY STARTING: " + n + " images");
        System.out.println("========================================");
        this.statusLabel.setText("Starting batch photometry on " + n + " images...");
        this.statusLabel.setForeground(Color.BLUE);
        SwingWorker<Void, String> swingWorker = new SwingWorker<Void, String>(this){
            final /* synthetic */ TargetPhotometryPanel this$0;
            {
                this.this$0 = targetPhotometryPanel;
            }

            /*
             * WARNING - Removed try catching itself - possible behaviour change.
             */
            @Override
            protected Void doInBackground() throws Exception {
                int n4 = 0;
                for (File file : arrayList) {
                    this.publish("Processing " + ++n4 + "/" + n + ": " + file.getName());
                    try {
                        Fits fits = new Fits(file);
                        ImageHDU imageHDU = (ImageHDU)fits.getHDU(0);
                        Header header = imageHDU.getHeader();
                        Object object = ((ImageData)imageHDU.getData()).getData();
                        int[][] nArray = this.this$0.convertToIntArray(object);
                        int n2 = Integer.MAX_VALUE;
                        int n3 = Integer.MIN_VALUE;
                        for (int i = 0; i < nArray.length; ++i) {
                            for (int j = 0; j < nArray[i].length; ++j) {
                                n2 = Math.min(n2, nArray[i][j]);
                                n3 = Math.max(n3, nArray[i][j]);
                            }
                        }
                        System.out.println("BATCH IMAGE: " + file.getName() + " - dimensions=" + nArray.length + "x" + nArray[0].length + ", min=" + n2 + ", max=" + n3);
                        fits.close();
                        List<PhotometryResult> list = this.this$0.processImagePhotometry(file, header, nArray);
                        TargetPhotometryPanel targetPhotometryPanel = this.this$0;
                        synchronized (targetPhotometryPanel) {
                            if (this.this$0.photometryResults == null) {
                                this.this$0.photometryResults = new ArrayList<PhotometryResult>();
                            }
                            this.this$0.photometryResults.addAll(list);
                        }
                        this.publish("Completed " + file.getName() + " (" + list.size() + " measurements)");
                    }
                    catch (Exception exception) {
                        this.publish("ERROR processing " + file.getName() + ": " + exception.getMessage());
                        exception.printStackTrace();
                    }
                }
                return null;
            }

            @Override
            protected void process(List<String> list) {
                for (String string : list) {
                    System.out.println("Batch Photometry: " + string);
                }
                if (!list.isEmpty()) {
                    this.this$0.statusLabel.setText(list.get(list.size() - 1));
                }
            }

            @Override
            protected void done() {
                try {
                    this.get();
                    this.this$0.updatePhotometryResultsTable();
                    int n2 = this.this$0.photometryResults != null ? this.this$0.photometryResults.size() : 0;
                    this.this$0.statusLabel.setText("Batch photometry complete - " + n2 + " total measurements from " + n + " images");
                    this.this$0.statusLabel.setForeground(new Color(0, 128, 0));
                    this.this$0.showApertures = true;
                    this.this$0.imagePanel.repaint();
                }
                catch (Exception exception) {
                    exception.printStackTrace();
                    this.this$0.statusLabel.setText("Batch photometry failed: " + exception.getMessage());
                    this.this$0.statusLabel.setForeground(Color.RED);
                }
            }
        };
        swingWorker.execute();
    }

    private List<PhotometryResult> processImagePhotometry(File file, Header header, int[][] nArray) throws Exception {
        Object object;
        Object object2;
        Object[] objectArray;
        double d;
        ArrayList<PhotometryResult> arrayList = new ArrayList<PhotometryResult>();
        String string = header.getStringValue("FILTER");
        if (string == null || string.trim().isEmpty()) {
            string = "Unknown";
        }
        int[][] nArray2 = this.currentImageData;
        this.currentImageData = nArray;
        this.selectCheckStar();
        this.currentImageData = nArray2;
        int n = this.currentRadialProfile.objectRadius;
        int n2 = this.currentRadialProfile.skyInnerRadius;
        int n3 = this.currentRadialProfile.skyOuterRadius;
        for (TargetInfo object3 : this.variableStarsInField) {
            int[] nArray3;
            if (object3.maxMag != null && object3.minMag != null) {
                try {
                    double nArray32 = Double.parseDouble(object3.maxMag.replaceAll("[^0-9.-]", ""));
                    double dArray = Double.parseDouble(object3.minMag.replaceAll("[^0-9.-]", ""));
                    d = Math.abs(nArray32 - dArray);
                    if (d < this.minVarDeltaV) {
                        continue;
                    }
                }
                catch (NumberFormatException string2) {
                    // empty catch block
                }
            }
            if ((nArray3 = this.raDecToPixel(object3.ra, object3.dec)) == null) {
                System.out.println("CENTROID DEBUG " + object3.name + " (VAR): raDecToPixel returned null");
                continue;
            }
            System.out.println("CENTROID DEBUG " + object3.name + " (VAR): catalog RA=" + String.format("%.6f", object3.ra) + ", Dec=" + String.format("%.6f", object3.dec) + " -> initial pixel (" + nArray3[0] + ", " + nArray3[1] + ")");
            objectArray = this.refineCentroid(nArray, nArray3[0], nArray3[1], 5);
            if (objectArray == null) {
                System.out.println("CENTROID DEBUG " + object3.name + " (VAR): refineCentroid failed");
                continue;
            }
            System.out.println("CENTROID DEBUG " + object3.name + " (VAR): refined centroid (" + String.format("%.2f", objectArray[0]) + ", " + String.format("%.2f", objectArray[1]) + ") shift=(" + String.format("%.2f", objectArray[0] - (double)nArray3[0]) + ", " + String.format("%.2f", objectArray[1] - (double)nArray3[1]) + ")");
            FluxMeasurement fluxMeasurement = this.measureFlux(nArray, objectArray[0], objectArray[1], n, n2, n3, object3.name);
            if (fluxMeasurement == null || fluxMeasurement.netFlux <= 0.0) continue;
            double d2 = -2.5 * Math.log10(fluxMeasurement.netFlux);
            double d3 = 1.0857 * fluxMeasurement.error / fluxMeasurement.netFlux;
            double[] dArray = this.pixelToRaDec(objectArray[0], objectArray[1]);
            object2 = file.getName();
            object = new PhotometryResult((String)object2, object3.name, "Variable", string, object3.ra, object3.dec, objectArray[0], objectArray[1], dArray != null ? dArray[0] : object3.ra, dArray != null ? dArray[1] : object3.dec, fluxMeasurement.starFlux, fluxMeasurement.skyFlux, fluxMeasurement.skyPerPixel, fluxMeasurement.netFlux, d2, d3);
            arrayList.add((PhotometryResult)object);
        }
        if (this.referenceStars != null) {
            for (ReferenceStar referenceStar : this.referenceStars) {
                Object object3;
                String string2 = referenceStar.label.isEmpty() ? referenceStar.name : referenceStar.label;
                objectArray = this.raDecToPixel(referenceStar.ra, referenceStar.dec);
                if (objectArray == null) {
                    System.out.println("CENTROID DEBUG " + string2 + ": raDecToPixel returned null for RA=" + String.format("%.6f", referenceStar.ra) + ", Dec=" + String.format("%.6f", referenceStar.dec));
                    continue;
                }
                System.out.println("CENTROID DEBUG " + string2 + ": catalog RA=" + String.format("%.6f", referenceStar.ra) + ", Dec=" + String.format("%.6f", referenceStar.dec) + " -> initial pixel (" + objectArray[0] + ", " + objectArray[1] + ")");
                double[] dArray = this.refineCentroid(nArray, objectArray[0], objectArray[1], 5);
                if (dArray == null) {
                    System.out.println("CENTROID DEBUG " + string2 + ": refineCentroid failed");
                    continue;
                }
                System.out.println("CENTROID DEBUG " + string2 + ": refined centroid (" + String.format("%.2f", dArray[0]) + ", " + String.format("%.2f", dArray[1]) + ") shift=(" + String.format("%.2f", dArray[0] - (double)objectArray[0]) + ", " + String.format("%.2f", dArray[1] - (double)objectArray[1]) + ")");
                FluxMeasurement fluxMeasurement = this.measureFlux(nArray, dArray[0], dArray[1], n, n2, n3, string2);
                if (fluxMeasurement == null || fluxMeasurement.netFlux <= 0.0) {
                    System.out.println("CENTROID DEBUG " + string2 + ": flux measurement failed or netFlux <= 0");
                    continue;
                }
                d = -2.5 * Math.log10(fluxMeasurement.netFlux);
                double d4 = 1.0857 * fluxMeasurement.error / fluxMeasurement.netFlux;
                object2 = this.pixelToRaDec(dArray[0], dArray[1]);
                if (object2 != null) {
                    System.out.println("CENTROID DEBUG " + string2 + ": centroid pixel (" + String.format("%.2f", dArray[0]) + ", " + String.format("%.2f", dArray[1]) + ") -> RA=" + String.format("%.6f", (double)object2[0]) + ", Dec=" + String.format("%.6f", (double)object2[1]) + " (offset from catalog: RA=" + String.format("%.6f", (double)(object2[0] - referenceStar.ra)) + ", Dec=" + String.format("%.6f", (double)(object2[1] - referenceStar.dec)) + ")");
                } else {
                    System.out.println("CENTROID DEBUG " + string2 + ": pixelToRaDec returned null");
                }
                object = file.getName();
                String string3 = "Comparison";
                if (this.checkStar != null) {
                    Object object4 = object3 = this.checkStar.label.isEmpty() ? this.checkStar.name : this.checkStar.label;
                    if (string2.equals(object3)) {
                        string3 = "Check";
                        System.out.println("CHECK STAR IDENTIFIED: " + string2 + " - setting type to 'Check'");
                    }
                }
                object3 = new PhotometryResult((String)object, string2, string3, string, referenceStar.ra, referenceStar.dec, dArray[0], dArray[1], object2 != null ? (double)object2[0] : referenceStar.ra, object2 != null ? (double)object2[1] : referenceStar.dec, fluxMeasurement.starFlux, fluxMeasurement.skyFlux, fluxMeasurement.skyPerPixel, fluxMeasurement.netFlux, d, d4);
                arrayList.add((PhotometryResult)object3);
            }
        }
        if (this.referenceStars != null && !this.referenceStars.isEmpty()) {
            File file2 = this.currentFitsFile;
            this.currentFitsFile = file;
            this.calculateZeropoints(arrayList);
            this.currentFitsFile = file2;
        }
        return arrayList;
    }

    private int[] raDecToPixel(double d, double d2) {
        try {
            Double d3 = this.currentFitsHeader.getDoubleValue("CRVAL1");
            Double d4 = this.currentFitsHeader.getDoubleValue("CRVAL2");
            Double d5 = this.currentFitsHeader.getDoubleValue("CRPIX1");
            Double d6 = this.currentFitsHeader.getDoubleValue("CRPIX2");
            if (d3 == null || d4 == null || d5 == null || d6 == null) {
                return null;
            }
            Double d7 = this.currentFitsHeader.getDoubleValue("CD1_1");
            Double d8 = this.currentFitsHeader.getDoubleValue("CD1_2");
            Double d9 = this.currentFitsHeader.getDoubleValue("CD2_1");
            Double d10 = this.currentFitsHeader.getDoubleValue("CD2_2");
            if (d7 == null || d8 == null || d9 == null || d10 == null) {
                Double d11 = this.currentFitsHeader.getDoubleValue("CDELT1");
                Double d12 = this.currentFitsHeader.getDoubleValue("CDELT2");
                Double d13 = this.currentFitsHeader.getDoubleValue("CROTA2");
                if (d11 != null && d12 != null) {
                    double d14 = d13 != null ? Math.toRadians(d13) : 0.0;
                    d7 = d11 * Math.cos(d14);
                    d8 = -d12.doubleValue() * Math.sin(d14);
                    d9 = d11 * Math.sin(d14);
                    d10 = d12 * Math.cos(d14);
                } else {
                    return null;
                }
            }
            double d15 = d - d3;
            double d16 = d2 - d4;
            d15 *= Math.cos(Math.toRadians(d4));
            double d17 = d7 * d10 - d8 * d9;
            if (Math.abs(d17) < 1.0E-10) {
                return null;
            }
            double d18 = (d10 * d15 - d8 * d16) / d17;
            double d19 = (-d9.doubleValue() * d15 + d7 * d16) / d17;
            int n = (int)Math.round(d5 + d18 - 1.0);
            int n2 = (int)Math.round(d6 + d19 - 1.0);
            return new int[]{n, n2};
        }
        catch (Exception exception) {
            return null;
        }
    }

    private double[] pixelToRaDec(double d, double d2) {
        try {
            Double d3 = this.currentFitsHeader.getDoubleValue("CRVAL1");
            Double d4 = this.currentFitsHeader.getDoubleValue("CRVAL2");
            Double d5 = this.currentFitsHeader.getDoubleValue("CRPIX1");
            Double d6 = this.currentFitsHeader.getDoubleValue("CRPIX2");
            if (d3 == null || d4 == null || d5 == null || d6 == null) {
                return null;
            }
            Double d7 = this.currentFitsHeader.getDoubleValue("CD1_1");
            Double d8 = this.currentFitsHeader.getDoubleValue("CD1_2");
            Double d9 = this.currentFitsHeader.getDoubleValue("CD2_1");
            Double d10 = this.currentFitsHeader.getDoubleValue("CD2_2");
            if (d7 == null || d8 == null || d9 == null || d10 == null) {
                Double d11 = this.currentFitsHeader.getDoubleValue("CDELT1");
                Double d12 = this.currentFitsHeader.getDoubleValue("CDELT2");
                Double d13 = this.currentFitsHeader.getDoubleValue("CROTA2");
                if (d11 != null && d12 != null) {
                    double d14 = d13 != null ? Math.toRadians(d13) : 0.0;
                    d7 = d11 * Math.cos(d14);
                    d8 = -d12.doubleValue() * Math.sin(d14);
                    d9 = d11 * Math.sin(d14);
                    d10 = d12 * Math.cos(d14);
                } else {
                    return null;
                }
            }
            double d15 = d - d5 + 1.0;
            double d16 = d2 - d6 + 1.0;
            double d17 = d7 * d15 + d8 * d16;
            double d18 = d9 * d15 + d10 * d16;
            double d19 = d3 + (d17 /= Math.cos(Math.toRadians(d4)));
            double d20 = d4 + d18;
            return new double[]{d19, d20};
        }
        catch (Exception exception) {
            return null;
        }
    }

    private double[] refineCentroid(int[][] nArray, int n, int n2, int n3) {
        try {
            int n4 = nArray.length;
            int n5 = nArray[0].length;
            int n6 = Math.max(0, n - n3);
            int n7 = Math.min(n4 - 1, n + n3);
            int n8 = Math.max(0, n2 - n3);
            int n9 = Math.min(n5 - 1, n2 + n3);
            int n10 = Integer.MAX_VALUE;
            for (int i = n6; i <= n7; ++i) {
                for (int j = n8; j <= n9; ++j) {
                    n10 = Math.min(n10, nArray[j][i]);
                }
            }
            double d = 0.0;
            double d2 = 0.0;
            double d3 = 0.0;
            for (int i = n6; i <= n7; ++i) {
                for (int j = n8; j <= n9; ++j) {
                    int n11 = nArray[j][i] - n10;
                    if (n11 <= 0) continue;
                    d += (double)(i * n11);
                    d2 += (double)(j * n11);
                    d3 += (double)n11;
                }
            }
            if (d3 == 0.0) {
                return null;
            }
            double d4 = d / d3;
            double d5 = d2 / d3;
            return new double[]{d4, d5};
        }
        catch (Exception exception) {
            return null;
        }
    }

    private FluxMeasurement measureFlux(int[][] nArray, double d, double d2, int n, int n2, int n3) {
        return this.measureFlux(nArray, d, d2, n, n2, n3, null);
    }

    private FluxMeasurement measureFlux(int[][] nArray, double d, double d2, int n, int n2, int n3, String string) {
        try {
            int n4;
            int n5 = nArray.length;
            int n6 = nArray[0].length;
            int n7 = (int)Math.round(d);
            int n8 = (int)Math.round(d2);
            if (n7 < n3 || n8 < n3 || n7 >= n5 - n3 || n8 >= n6 - n3) {
                if (string != null) {
                    System.out.println("DEBUG " + string + ": flux failed - bounds check (cx=" + n7 + ", cy=" + n8 + ")");
                }
                return null;
            }
            double d3 = 0.0;
            int n9 = 0;
            int n10 = 0;
            for (int i = n7 - n; i <= n7 + n; ++i) {
                for (n4 = n8 - n; n4 <= n8 + n; ++n4) {
                    double d4 = Math.sqrt(((double)i - d) * ((double)i - d) + ((double)n4 - d2) * ((double)n4 - d2));
                    if (!(d4 <= (double)n)) continue;
                    int n11 = nArray[n4][i];
                    d3 += (double)n11;
                    ++n9;
                    if (n11 <= n10) continue;
                    n10 = n11;
                }
            }
            if (string != null && (string.equals("123") || string.equals("129"))) {
                System.out.println("DEBUG " + string + ": star aperture max pixel value = " + n10);
                System.out.println("DEBUG " + string + ": star aperture average = " + d3 / (double)n9);
            }
            ArrayList<Integer> arrayList = new ArrayList<Integer>();
            for (n4 = n7 - n3; n4 <= n7 + n3; ++n4) {
                for (int i = n8 - n3; i <= n8 + n3; ++i) {
                    double d5 = Math.sqrt(((double)n4 - d) * ((double)n4 - d) + ((double)i - d2) * ((double)i - d2));
                    if (!(d5 >= (double)n2) || !(d5 <= (double)n3)) continue;
                    arrayList.add(nArray[i][n4]);
                }
            }
            if (arrayList.isEmpty()) {
                if (string != null) {
                    System.out.println("DEBUG " + string + ": flux failed - no sky pixels");
                }
                return null;
            }
            Collections.sort(arrayList);
            n4 = arrayList.size();
            double d6 = n4 % 2 == 0 ? (double)((Integer)arrayList.get(n4 / 2 - 1) + (Integer)arrayList.get(n4 / 2)) / 2.0 : (double)((Integer)arrayList.get(n4 / 2)).intValue();
            double d7 = 0.0;
            Iterator iterator = arrayList.iterator();
            while (iterator.hasNext()) {
                int n12 = (Integer)iterator.next();
                d7 += (double)n12;
            }
            double d8 = d7 / (double)n4;
            double d9 = 3.0 * d6 - 2.0 * d8;
            if (d9 < 0.0) {
                d9 = d6;
            }
            if (d9 > d6) {
                d9 = d6;
            }
            if (string != null && (string.equals("123") || string.equals("129"))) {
                System.out.println("DEBUG " + string + ": sky annulus has " + n4 + " pixels");
                System.out.println("DEBUG " + string + ": sky min=" + String.valueOf(arrayList.get(0)) + ", max=" + String.valueOf(arrayList.get(n4 - 1)) + ", median=" + d6 + ", mean=" + d8);
                System.out.println("DEBUG " + string + ": sky Q1=" + String.valueOf(arrayList.get(n4 / 4)) + ", Q3=" + String.valueOf(arrayList.get(3 * n4 / 4)));
                System.out.println("DEBUG " + string + ": skyPerPixel (mode) = " + d9);
            }
            double d10 = d9 * (double)n4;
            double d11 = d9 * (double)n9;
            double d12 = d3 - d11;
            if (string != null) {
                System.out.println("DEBUG " + string + ": starSum=" + d3 + ", skyPerPixel=" + d9 + ", starPixels=" + n9 + ", netFlux=" + d12);
            }
            if (d12 <= 0.0) {
                if (string != null) {
                    System.out.println("DEBUG " + string + ": flux failed - netFlux <= 0 (" + d12 + ")");
                }
                return null;
            }
            double d13 = Math.sqrt(d12 + (double)n9 * d9);
            return new FluxMeasurement(d3, d10, d9, d12, d13);
        }
        catch (Exception exception) {
            return null;
        }
    }

    private void updatePhotometryResultsTable() {
        Object object;
        if (this.photometryResults == null || this.photometryResults.isEmpty()) {
            return;
        }
        this.photometryResultsPanel.removeAll();
        this.photometryResultsPanel.setLayout(new BorderLayout(10, 10));
        String string = !this.photometryResults.isEmpty() && this.photometryResults.get((int)0).filter != null ? this.photometryResults.get((int)0).filter : "Std";
        Object[] objectArray = new String[]{"Image ID", "Name", "Type", "Cat RA", "Cat Dec", "Centroid X", "Centroid Y", "Centroid RA", "Centroid Dec", "Net Flux", "Inst Mag", "Inst Err", "Std " + string, "Std Err"};
        Object[][] objectArray2 = new Object[this.photometryResults.size()][objectArray.length];
        for (int i = 0; i < this.photometryResults.size(); ++i) {
            object = this.photometryResults.get(i);
            objectArray2[i][0] = ((PhotometryResult)object).imageId;
            objectArray2[i][1] = ((PhotometryResult)object).name;
            objectArray2[i][2] = ((PhotometryResult)object).type;
            objectArray2[i][3] = String.format("%.6f", ((PhotometryResult)object).catalogRA);
            objectArray2[i][4] = String.format("%.6f", ((PhotometryResult)object).catalogDec);
            objectArray2[i][5] = String.format("%.2f", ((PhotometryResult)object).centroidX);
            objectArray2[i][6] = String.format("%.2f", ((PhotometryResult)object).centroidY);
            objectArray2[i][7] = String.format("%.6f", ((PhotometryResult)object).centroidRA);
            objectArray2[i][8] = String.format("%.6f", ((PhotometryResult)object).centroidDec);
            objectArray2[i][9] = String.format("%.1f", ((PhotometryResult)object).netFlux);
            objectArray2[i][10] = String.format("%.3f", ((PhotometryResult)object).instMag);
            objectArray2[i][11] = String.format("%.3f", ((PhotometryResult)object).magError);
            if (!Double.isNaN(((PhotometryResult)object).calibratedMag)) {
                objectArray2[i][12] = String.format("%.3f", ((PhotometryResult)object).calibratedMag);
                objectArray2[i][13] = String.format("%.3f", ((PhotometryResult)object).magError);
                continue;
            }
            objectArray2[i][12] = "--";
            objectArray2[i][13] = "--";
        }
        JTable jTable = new JTable(objectArray2, objectArray);
        jTable.setAutoResizeMode(0);
        jTable.setFont(new Font("Monospaced", 0, 11));
        jTable.getColumnModel().getColumn(0).setPreferredWidth(200);
        jTable.getColumnModel().getColumn(1).setPreferredWidth(120);
        jTable.getColumnModel().getColumn(2).setPreferredWidth(90);
        jTable.getColumnModel().getColumn(3).setPreferredWidth(100);
        jTable.getColumnModel().getColumn(4).setPreferredWidth(100);
        jTable.getColumnModel().getColumn(5).setPreferredWidth(90);
        jTable.getColumnModel().getColumn(6).setPreferredWidth(90);
        jTable.getColumnModel().getColumn(7).setPreferredWidth(100);
        jTable.getColumnModel().getColumn(8).setPreferredWidth(100);
        jTable.getColumnModel().getColumn(9).setPreferredWidth(100);
        jTable.getColumnModel().getColumn(10).setPreferredWidth(80);
        jTable.getColumnModel().getColumn(11).setPreferredWidth(80);
        jTable.getColumnModel().getColumn(12).setPreferredWidth(80);
        jTable.getColumnModel().getColumn(13).setPreferredWidth(80);
        object = new JScrollPane(jTable);
        this.photometryResultsPanel.add((Component)object, "Center");
        int n = (int)this.photometryResults.stream().filter(photometryResult -> photometryResult.type.equals("Variable")).count();
        int n2 = (int)this.photometryResults.stream().filter(photometryResult -> photometryResult.type.equals("Comparison")).count();
        int n3 = (int)this.photometryResults.stream().filter(photometryResult -> photometryResult.type.equals("Check")).count();
        JLabel jLabel = new JLabel(String.format("Total: %d stars (%d variables, %d comparisons, %d check)", this.photometryResults.size(), n, n2, n3), 0);
        jLabel.setFont(new Font("SansSerif", 1, 12));
        jLabel.setBorder(BorderFactory.createEmptyBorder(10, 10, 10, 10));
        this.photometryResultsPanel.add((Component)jLabel, "South");
        this.photometryResultsPanel.revalidate();
        this.photometryResultsPanel.repaint();
        this.updateVariableSelector();
    }

    private void handleStarClick(int n, int n2) {
        int n3;
        int n4;
        int[] nArray;
        if (this.currentImage == null || this.currentFitsHeader == null) {
            return;
        }
        this.selectedStar = null;
        if (this.variableStarsInField != null) {
            for (TargetInfo object : this.variableStarsInField) {
                nArray = this.getStarScreenCoords(object.ra, object.dec);
                if (nArray == null || !(Math.sqrt((n4 = n - nArray[0]) * n4 + (n3 = n2 - nArray[1]) * n3) < 10.0)) continue;
                this.selectedStar = object;
                this.selectedStarRA = object.ra;
                this.selectedStarDec = object.dec;
                this.showStarDetails(object);
                this.imagePanel.repaint();
                return;
            }
        }
        if (this.referenceStars != null) {
            for (ReferenceStar referenceStar : this.referenceStars) {
                nArray = this.getStarScreenCoords(referenceStar.ra, referenceStar.dec);
                if (nArray == null || !(Math.sqrt((n4 = n - nArray[0]) * n4 + (n3 = n2 - nArray[1]) * n3) < 10.0)) continue;
                this.selectedStar = referenceStar;
                this.selectedStarRA = referenceStar.ra;
                this.selectedStarDec = referenceStar.dec;
                this.showStarDetails(referenceStar);
                this.imagePanel.repaint();
                return;
            }
        }
        this.selectedStar = null;
        this.detailsArea.setText("");
        this.imagePanel.repaint();
    }

    private int[] getStarScreenCoords(double d, double d2) {
        try {
            Double d3 = this.currentFitsHeader.getDoubleValue("CRVAL1");
            Double d4 = this.currentFitsHeader.getDoubleValue("CRVAL2");
            Double d5 = this.currentFitsHeader.getDoubleValue("CRPIX1");
            Double d6 = this.currentFitsHeader.getDoubleValue("CRPIX2");
            if (d3 == null || d4 == null || d5 == null || d6 == null) {
                return null;
            }
            Double d7 = this.currentFitsHeader.getDoubleValue("CD1_1");
            Double d8 = this.currentFitsHeader.getDoubleValue("CD1_2");
            Double d9 = this.currentFitsHeader.getDoubleValue("CD2_1");
            Double d10 = this.currentFitsHeader.getDoubleValue("CD2_2");
            if (d7 == null || d8 == null || d9 == null || d10 == null) {
                Double d11 = this.currentFitsHeader.getDoubleValue("CDELT1");
                Double d12 = this.currentFitsHeader.getDoubleValue("CDELT2");
                Double d13 = this.currentFitsHeader.getDoubleValue("CROTA2");
                if (d11 != null && d12 != null) {
                    double d14 = d13 != null ? Math.toRadians(d13) : 0.0;
                    d7 = d11 * Math.cos(d14);
                    d8 = -d12.doubleValue() * Math.sin(d14);
                    d9 = d11 * Math.sin(d14);
                    d10 = d12 * Math.cos(d14);
                } else {
                    return null;
                }
            }
            double d15 = d - d3;
            double d16 = d2 - d4;
            d15 *= Math.cos(Math.toRadians(d4));
            double d17 = d7 * d10 - d8 * d9;
            if (Math.abs(d17) < 1.0E-10) {
                return null;
            }
            double d18 = (d10 * d15 - d8 * d16) / d17;
            double d19 = (-d9.doubleValue() * d15 + d7 * d16) / d17;
            double d20 = d5 + d18 - 1.0;
            double d21 = d6 + d19 - 1.0;
            int n = this.imagePanel.getWidth();
            int n2 = this.imagePanel.getHeight();
            int n3 = this.currentImage.getWidth();
            int n4 = this.currentImage.getHeight();
            double d22 = (double)n / (double)n3;
            double d23 = (double)n2 / (double)n4;
            double d24 = Math.min(d22, d23) * this.zoomLevel;
            int n5 = (int)((double)n3 * d24);
            int n6 = (int)((double)n4 * d24);
            int n7 = (n - n5) / 2 + this.panOffsetX;
            int n8 = (n2 - n6) / 2 + this.panOffsetY;
            int n9 = n7 + (int)(d20 * d24);
            int n10 = n8 + (int)(d21 * d24);
            return new int[]{n9, n10};
        }
        catch (Exception exception) {
            return null;
        }
    }

    private String raToSexagesimal(double d) {
        int n = (int)(d / 15.0);
        double d2 = (d / 15.0 - (double)n) * 60.0;
        int n2 = (int)d2;
        double d3 = (d2 - (double)n2) * 60.0;
        return String.format("%02dh%02dm%04.1fs", n, n2, d3);
    }

    private String decToSexagesimal(double d) {
        String string = d >= 0.0 ? "+" : "-";
        d = Math.abs(d);
        int n = (int)d;
        double d2 = (d - (double)n) * 60.0;
        int n2 = (int)d2;
        double d3 = (d2 - (double)n2) * 60.0;
        return String.format("%s%02d\u00b0%02d'%04.1f\"", string, n, n2, d3);
    }

    private void showStarDetails(Object object) {
        StringBuilder stringBuilder = new StringBuilder();
        if (object instanceof TargetInfo) {
            TargetInfo targetInfo = (TargetInfo)object;
            stringBuilder.append(targetInfo.name).append(" | ");
            stringBuilder.append("RA: ").append(this.raToSexagesimal(targetInfo.ra)).append(" | ");
            stringBuilder.append("Dec: ").append(this.decToSexagesimal(targetInfo.dec)).append(" | ");
            if (targetInfo.maxMag != null && targetInfo.minMag != null) {
                stringBuilder.append("Mag: ").append(targetInfo.minMag).append("-").append(targetInfo.maxMag).append(" | ");
            }
            if (targetInfo.period > 0.0) {
                stringBuilder.append(String.format("Period: %.4f d", targetInfo.period)).append(" | ");
            }
            if (targetInfo.variabilityType != null && !targetInfo.variabilityType.isEmpty()) {
                stringBuilder.append("Type: ").append(targetInfo.variabilityType).append(" | ");
            }
            if (targetInfo.auid != null && !targetInfo.auid.isEmpty()) {
                stringBuilder.append("AUID: ").append(targetInfo.auid);
            }
        } else if (object instanceof ReferenceStar) {
            ReferenceStar referenceStar = (ReferenceStar)object;
            stringBuilder.append("Comp ").append(referenceStar.name).append(" | ");
            stringBuilder.append("RA: ").append(this.raToSexagesimal(referenceStar.ra)).append(" | ");
            stringBuilder.append("Dec: ").append(this.decToSexagesimal(referenceStar.dec)).append(" | ");
            stringBuilder.append(String.format("Mag V: %.2f", referenceStar.mag));
            if (referenceStar.magError > 0.0) {
                stringBuilder.append(String.format("(%.2f)", referenceStar.magError));
            }
            stringBuilder.append(" | ");
            if (referenceStar.bv != 0.0) {
                stringBuilder.append(String.format("B-V: %.2f", referenceStar.bv));
                if (referenceStar.bvError > 0.0) {
                    stringBuilder.append(String.format("(%.2f)", referenceStar.bvError));
                }
            }
        }
        this.detailsArea.setText(stringBuilder.toString());
    }

    private void updateVariableSelector() {
        if (this.variableSelector == null || this.photometryResults == null) {
            return;
        }
        HashSet<String> hashSet = new HashSet<String>();
        for (PhotometryResult object2 : this.photometryResults) {
            if (!object2.type.equals("Variable") || Double.isNaN(object2.calibratedMag)) continue;
            hashSet.add(object2.name);
        }
        String string = (String)this.variableSelector.getSelectedItem();
        this.variableSelector.removeAllItems();
        ArrayList arrayList = new ArrayList(hashSet);
        Collections.sort(arrayList);
        for (String string2 : arrayList) {
            this.variableSelector.addItem(string2);
        }
        if (string != null && hashSet.contains(string)) {
            this.variableSelector.setSelectedItem(string);
        } else if (!arrayList.isEmpty()) {
            this.variableSelector.setSelectedIndex(0);
            this.updateLightcurve();
        }
    }

    private static class StarLocation {
        int x;
        int y;
        int peakValue;

        StarLocation(int n, int n2, int n3) {
            this.x = n;
            this.y = n2;
            this.peakValue = n3;
        }
    }

    private static class RadialProfile {
        double[] intensity;
        int objectRadius;
        int skyInnerRadius;
        int skyOuterRadius;
        int centerX;
        int centerY;

        RadialProfile(double[] dArray, int n, int n2, int n3, int n4, int n5) {
            this.intensity = dArray;
            this.objectRadius = n;
            this.skyInnerRadius = n2;
            this.skyOuterRadius = n3;
            this.centerX = n4;
            this.centerY = n5;
        }
    }

    private static class ReferenceStar {
        String name;
        double ra;
        double dec;
        double mag;
        double magError;
        double bv;
        double bvError;
        String label;

        ReferenceStar(String string, double d, double d2, double d3) {
            this.name = string;
            this.ra = d;
            this.dec = d2;
            this.mag = d3;
            this.magError = 0.0;
            this.bv = 0.0;
            this.bvError = 0.0;
            this.label = "";
        }

        ReferenceStar(String string, double d, double d2, double d3, double d4, double d5, double d6, String string2) {
            this.name = string;
            this.ra = d;
            this.dec = d2;
            this.mag = d3;
            this.magError = d4;
            this.bv = d5;
            this.bvError = d6;
            this.label = string2 != null ? string2 : "";
        }
    }

    private static class StatusCellRenderer
    extends DefaultTableCellRenderer {
        private StatusCellRenderer() {
        }

        @Override
        public Component getTableCellRendererComponent(JTable jTable, Object object, boolean bl, boolean bl2, int n, int n2) {
            Component component = super.getTableCellRendererComponent(jTable, object, bl, bl2, n, n2);
            if (!bl && object != null) {
                String string;
                switch (string = object.toString()) {
                    case "Not Done": {
                        component.setBackground(new Color(255, 200, 200));
                        break;
                    }
                    case "Done": {
                        component.setBackground(new Color(200, 255, 200));
                        break;
                    }
                    case "Bad Obs": {
                        component.setBackground(new Color(200, 200, 200));
                        break;
                    }
                    default: {
                        component.setBackground(Color.WHITE);
                    }
                }
            }
            return component;
        }
    }

    private class ImageDisplayPanel
    extends JPanel {
        private BufferedImage image;
        private Header fitsHeader;
        private List<ReferenceStar> stars;
        private List<TargetInfo> variableStars;
        private String message;
        private double northAngle = 0.0;
        private double eastAngle = 0.0;
        private int mouseX = -1;
        private int mouseY = -1;
        private int imageOffsetX = 0;
        private int imageOffsetY = 0;
        private double imageScale = 1.0;

        public ImageDisplayPanel() {
            this.setPreferredSize(new Dimension(800, 800));
            this.setBackground(Color.BLACK);
            this.addMouseMotionListener(new MouseMotionAdapter(){

                @Override
                public void mouseMoved(MouseEvent mouseEvent) {
                    ImageDisplayPanel.this.mouseX = mouseEvent.getX();
                    ImageDisplayPanel.this.mouseY = mouseEvent.getY();
                    ImageDisplayPanel.this.updateMouseInfo(ImageDisplayPanel.this.mouseX, ImageDisplayPanel.this.mouseY);
                    ImageDisplayPanel.this.repaint();
                }

                @Override
                public void mouseDragged(MouseEvent mouseEvent) {
                    if (TargetPhotometryPanel.this.isDragging) {
                        int n = mouseEvent.getX() - TargetPhotometryPanel.this.dragStartX;
                        int n2 = mouseEvent.getY() - TargetPhotometryPanel.this.dragStartY;
                        TargetPhotometryPanel.this.panOffsetX += n;
                        TargetPhotometryPanel.this.panOffsetY += n2;
                        TargetPhotometryPanel.this.dragStartX = mouseEvent.getX();
                        TargetPhotometryPanel.this.dragStartY = mouseEvent.getY();
                        ImageDisplayPanel.this.updateMouseInfo(mouseEvent.getX(), mouseEvent.getY());
                        ImageDisplayPanel.this.repaint();
                    }
                }
            });
            this.addMouseListener(new MouseAdapter(){

                @Override
                public void mouseExited(MouseEvent mouseEvent) {
                    ImageDisplayPanel.this.mouseX = -1;
                    ImageDisplayPanel.this.mouseY = -1;
                    if (TargetPhotometryPanel.this.mouseInfoLabel != null) {
                        TargetPhotometryPanel.this.mouseInfoLabel.setText("x: -, y: -, intensity: -");
                    }
                    ImageDisplayPanel.this.repaint();
                }

                @Override
                public void mousePressed(MouseEvent mouseEvent) {
                    TargetPhotometryPanel.this.dragStartX = mouseEvent.getX();
                    TargetPhotometryPanel.this.dragStartY = mouseEvent.getY();
                    TargetPhotometryPanel.this.isDragging = true;
                }

                @Override
                public void mouseReleased(MouseEvent mouseEvent) {
                    TargetPhotometryPanel.this.isDragging = false;
                }

                @Override
                public void mouseClicked(MouseEvent mouseEvent) {
                    TargetPhotometryPanel.this.handleStarClick(mouseEvent.getX(), mouseEvent.getY());
                }
            });
            this.addMouseWheelListener(mouseWheelEvent -> {
                int n = mouseWheelEvent.getX();
                int n2 = mouseWheelEvent.getY();
                int n3 = this.getWidth();
                int n4 = this.getHeight();
                if (this.image != null) {
                    int n5 = this.image.getWidth();
                    int n6 = this.image.getHeight();
                    double d = (double)n3 / (double)n5;
                    double d2 = (double)n4 / (double)n6;
                    double d3 = Math.min(d, d2);
                    double d4 = d3 * TargetPhotometryPanel.this.zoomLevel;
                    int n7 = (int)((double)n5 * d4);
                    int n8 = (int)((double)n6 * d4);
                    int n9 = (n3 - n7) / 2 + TargetPhotometryPanel.this.panOffsetX;
                    int n10 = (n4 - n8) / 2 + TargetPhotometryPanel.this.panOffsetY;
                    double d5 = (double)(n - n9) / d4;
                    double d6 = (double)(n2 - n10) / d4;
                    TargetPhotometryPanel.this.zoomLevel = mouseWheelEvent.getWheelRotation() < 0 ? (TargetPhotometryPanel.this.zoomLevel *= 1.1) : (TargetPhotometryPanel.this.zoomLevel /= 1.1);
                    TargetPhotometryPanel.this.zoomLevel = Math.max(0.1, Math.min(10.0, TargetPhotometryPanel.this.zoomLevel));
                    double d7 = d3 * TargetPhotometryPanel.this.zoomLevel;
                    int n11 = (int)((double)n5 * d7);
                    int n12 = (int)((double)n6 * d7);
                    int n13 = (n3 - n11) / 2;
                    int n14 = (n4 - n12) / 2;
                    int n15 = n13 + (int)(d5 * d7);
                    int n16 = n14 + (int)(d6 * d7);
                    TargetPhotometryPanel.this.panOffsetX = n - n15;
                    TargetPhotometryPanel.this.panOffsetY = n2 - n16;
                }
                this.repaint();
            });
        }

        public void setImage(BufferedImage bufferedImage, Header header, List<ReferenceStar> list) {
            this.image = bufferedImage;
            this.fitsHeader = header;
            this.stars = list;
            this.variableStars = TargetPhotometryPanel.this.variableStarsInField;
            this.message = null;
            if (header != null) {
                this.calculateOrientation(header);
            }
            this.repaint();
        }

        public void setMessage(String string) {
            this.message = string;
            this.repaint();
        }

        private void calculateOrientation(Header header) {
            try {
                Double d = header.getDoubleValue("CD1_1");
                Double d2 = header.getDoubleValue("CD1_2");
                Double d3 = header.getDoubleValue("CD2_1");
                Double d4 = header.getDoubleValue("CD2_2");
                if (d != null && d2 != null && d3 != null && d4 != null) {
                    this.northAngle = Math.atan2(d2, d4);
                    this.eastAngle = this.northAngle - 1.5707963267948966;
                    return;
                }
                Double d5 = header.getDoubleValue("CROTA2");
                if (d5 != null) {
                    this.northAngle = Math.toRadians(d5);
                    this.eastAngle = this.northAngle - 1.5707963267948966;
                    return;
                }
                this.northAngle = 0.0;
                this.eastAngle = -1.5707963267948966;
            }
            catch (Exception exception) {
                this.northAngle = 0.0;
                this.eastAngle = -1.5707963267948966;
            }
        }

        @Override
        protected void paintComponent(Graphics graphics) {
            Object object;
            super.paintComponent(graphics);
            Graphics2D graphics2D = (Graphics2D)graphics;
            graphics2D.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
            int n = this.getWidth();
            int n2 = this.getHeight();
            if (this.message != null) {
                graphics2D.setColor(Color.WHITE);
                graphics2D.setFont(new Font("SansSerif", 0, 14));
                FontMetrics fontMetrics = graphics2D.getFontMetrics();
                int n3 = fontMetrics.stringWidth(this.message);
                graphics2D.drawString(this.message, (n - n3) / 2, n2 / 2);
                return;
            }
            if (this.image == null) {
                return;
            }
            int n4 = this.image.getWidth();
            int n5 = this.image.getHeight();
            double d = (double)n / (double)n4;
            double d2 = (double)n2 / (double)n5;
            double d3 = Math.min(d, d2) * TargetPhotometryPanel.this.zoomLevel;
            int n6 = (int)((double)n4 * d3);
            int n7 = (int)((double)n5 * d3);
            int n8 = (n - n6) / 2 + TargetPhotometryPanel.this.panOffsetX;
            int n9 = (n2 - n7) / 2 + TargetPhotometryPanel.this.panOffsetY;
            this.imageOffsetX = n8;
            this.imageOffsetY = n9;
            this.imageScale = d3;
            graphics2D.drawImage(this.image, n8, n9, n6, n7, null);
            this.drawCompass(graphics2D, n8 + 50, n9 + 50, 30, d3);
            if (this.variableStars != null && this.fitsHeader != null) {
                this.drawVariableStars(graphics2D, n8, n9, d3);
            }
            if (this.stars != null && this.fitsHeader != null) {
                this.drawReferenceStars(graphics2D, n8, n9, d3);
            }
            if (TargetPhotometryPanel.this.selectedStar != null && this.fitsHeader != null && (object = this.getStarScreenCoordsFromHeader(TargetPhotometryPanel.this.selectedStarRA, TargetPhotometryPanel.this.selectedStarDec, this.fitsHeader, n8, n9, d3)) != null) {
                graphics2D.setColor(Color.RED);
                float[] fArray = new float[]{5.0f, 5.0f};
                graphics2D.setStroke(new BasicStroke(1.0f, 0, 0, 10.0f, fArray, 0.0f));
                graphics2D.drawLine(0, object[1], n, object[1]);
                graphics2D.drawLine(object[0], 0, object[0], n2);
            }
            if (TargetPhotometryPanel.this.radialProfileStar != null) {
                int n10 = n8 + (int)((double)TargetPhotometryPanel.this.radialProfileStar.x * d3);
                int n11 = n9 + (int)((double)TargetPhotometryPanel.this.radialProfileStar.y * d3);
                graphics2D.setColor(new Color(138, 43, 226));
                graphics2D.setStroke(new BasicStroke(3.0f));
                int n12 = n10 - 20;
                int n13 = n11 - 20;
                int n14 = n10 - 8;
                int n15 = n11 - 8;
                graphics2D.drawLine(n12, n13, n14, n15);
                int n16 = 10;
                double d4 = Math.atan2(n15 - n13, n14 - n12);
                int[] nArray = new int[]{n14, n14 - (int)((double)n16 * Math.cos(d4 - 0.5235987755982988)), n14 - (int)((double)n16 * Math.cos(d4 + 0.5235987755982988))};
                int[] nArray2 = new int[]{n15, n15 - (int)((double)n16 * Math.sin(d4 - 0.5235987755982988)), n15 - (int)((double)n16 * Math.sin(d4 + 0.5235987755982988))};
                graphics2D.fillPolygon(nArray, nArray2, 3);
            }
            if (this.mouseX >= 0 && this.mouseY >= 0 && this.fitsHeader != null) {
                this.drawCoordinateBox(graphics2D, this.mouseX, this.mouseY);
            }
            if (TargetPhotometryPanel.this.currentFitsFile != null) {
                graphics2D.setFont(new Font("SansSerif", 1, 12));
                graphics2D.setColor(Color.YELLOW);
                object = TargetPhotometryPanel.this.currentFitsFile.getName();
                graphics2D.drawString((String)object, 10, 20);
            }
        }

        private int[] getStarScreenCoordsFromHeader(double d, double d2, Header header, int n, int n2, double d3) {
            try {
                Double d4 = header.getDoubleValue("CRVAL1");
                Double d5 = header.getDoubleValue("CRVAL2");
                Double d6 = header.getDoubleValue("CRPIX1");
                Double d7 = header.getDoubleValue("CRPIX2");
                if (d4 == null || d5 == null || d6 == null || d7 == null) {
                    return null;
                }
                Double d8 = header.getDoubleValue("CD1_1");
                Double d9 = header.getDoubleValue("CD1_2");
                Double d10 = header.getDoubleValue("CD2_1");
                Double d11 = header.getDoubleValue("CD2_2");
                if (d8 == null || d9 == null || d10 == null || d11 == null) {
                    Double d12 = header.getDoubleValue("CDELT1");
                    Double d13 = header.getDoubleValue("CDELT2");
                    Double d14 = header.getDoubleValue("CROTA2");
                    if (d12 != null && d13 != null) {
                        double d15 = d14 != null ? Math.toRadians(d14) : 0.0;
                        d8 = d12 * Math.cos(d15);
                        d9 = -d13.doubleValue() * Math.sin(d15);
                        d10 = d12 * Math.sin(d15);
                        d11 = d13 * Math.cos(d15);
                    } else {
                        return null;
                    }
                }
                double d16 = d - d4;
                double d17 = d2 - d5;
                d16 *= Math.cos(Math.toRadians(d5));
                double d18 = d8 * d11 - d9 * d10;
                if (Math.abs(d18) < 1.0E-10) {
                    return null;
                }
                double d19 = (d11 * d16 - d9 * d17) / d18;
                double d20 = (-d10.doubleValue() * d16 + d8 * d17) / d18;
                double d21 = d6 + d19 - 1.0;
                double d22 = d7 + d20 - 1.0;
                int n3 = n + (int)(d21 * d3);
                int n4 = n2 + (int)(d22 * d3);
                return new int[]{n3, n4};
            }
            catch (Exception exception) {
                return null;
            }
        }

        private void drawCompass(Graphics2D graphics2D, int n, int n2, int n3, double d) {
            graphics2D.setColor(Color.YELLOW);
            graphics2D.setStroke(new BasicStroke(2.0f));
            int n4 = n + (int)((double)n3 * Math.sin(this.northAngle));
            int n5 = n2 - (int)((double)n3 * Math.cos(this.northAngle));
            graphics2D.drawLine(n, n2, n4, n5);
            this.drawArrowHead(graphics2D, n, n2, n4, n5);
            graphics2D.setFont(new Font("SansSerif", 1, 14));
            graphics2D.drawString("N", n4 + 5, n5 + 5);
            int n6 = n + (int)((double)n3 * Math.sin(this.eastAngle));
            int n7 = n2 - (int)((double)n3 * Math.cos(this.eastAngle));
            graphics2D.drawLine(n, n2, n6, n7);
            this.drawArrowHead(graphics2D, n, n2, n6, n7);
            graphics2D.drawString("E", n6 + 5, n7 + 5);
        }

        private void drawArrowHead(Graphics2D graphics2D, int n, int n2, int n3, int n4) {
            double d = Math.atan2(n4 - n2, n3 - n);
            int n5 = 8;
            int[] nArray = new int[3];
            int[] nArray2 = new int[3];
            nArray[0] = n3;
            nArray2[0] = n4;
            nArray[1] = n3 - (int)((double)n5 * Math.cos(d - 0.5235987755982988));
            nArray2[1] = n4 - (int)((double)n5 * Math.sin(d - 0.5235987755982988));
            nArray[2] = n3 - (int)((double)n5 * Math.cos(d + 0.5235987755982988));
            nArray2[2] = n4 - (int)((double)n5 * Math.sin(d + 0.5235987755982988));
            graphics2D.fillPolygon(nArray, nArray2, 3);
        }

        private void drawVariableStars(Graphics2D graphics2D, int n, int n2, double d) {
            if (this.fitsHeader == null || this.variableStars == null || this.variableStars.isEmpty()) {
                return;
            }
            try {
                Double d2 = this.fitsHeader.getDoubleValue("CRVAL1");
                Double d3 = this.fitsHeader.getDoubleValue("CRVAL2");
                Double d4 = this.fitsHeader.getDoubleValue("CRPIX1");
                Double d5 = this.fitsHeader.getDoubleValue("CRPIX2");
                Double d6 = this.fitsHeader.getDoubleValue("CD1_1");
                Double d7 = this.fitsHeader.getDoubleValue("CD1_2");
                Double d8 = this.fitsHeader.getDoubleValue("CD2_1");
                Double d9 = this.fitsHeader.getDoubleValue("CD2_2");
                if (d2 == null || d3 == null || d4 == null || d5 == null) {
                    return;
                }
                if (d6 == null || d7 == null || d8 == null || d9 == null) {
                    Double d10 = this.fitsHeader.getDoubleValue("CDELT1");
                    Double object = this.fitsHeader.getDoubleValue("CDELT2");
                    Double exception = this.fitsHeader.getDoubleValue("CROTA2");
                    if (d10 != null && object != null) {
                        double d11 = exception != null ? Math.toRadians(exception) : 0.0;
                        d6 = d10 * Math.cos(d11);
                        d7 = -object.doubleValue() * Math.sin(d11);
                        d8 = d10 * Math.sin(d11);
                        d9 = object * Math.cos(d11);
                    } else {
                        return;
                    }
                }
                graphics2D.setColor(Color.YELLOW);
                graphics2D.setStroke(new BasicStroke(2.0f));
                for (TargetInfo targetInfo : this.variableStars) {
                    int n3;
                    double d12;
                    double d13;
                    if (targetInfo.maxMag != null && targetInfo.minMag != null) {
                        try {
                            double d16 = Double.parseDouble(targetInfo.maxMag.replaceAll("[^0-9.-]", ""));
                            d13 = Double.parseDouble(targetInfo.minMag.replaceAll("[^0-9.-]", ""));
                            d12 = Math.abs(d16 - d13);
                            if (d12 < TargetPhotometryPanel.this.minVarDeltaV) {
                                continue;
                            }
                        }
                        catch (Exception exception) {
                            // empty catch block
                        }
                    }
                    double d14 = targetInfo.ra - d2;
                    d13 = targetInfo.dec - d3;
                    d14 *= Math.cos(Math.toRadians(d3));
                    d12 = d6 * d9 - d7 * d8;
                    if (Math.abs(d12) < 1.0E-10) continue;
                    double d15 = (d9 * d14 - d7 * d13) / d12;
                    double d16 = (-d8.doubleValue() * d14 + d6 * d13) / d12;
                    double d17 = d4 + d15 - 1.0;
                    double d18 = d5 + d16 - 1.0;
                    int n4 = n + (int)(d17 * d);
                    int n5 = n2 + (int)(d18 * d);
                    if (TargetPhotometryPanel.this.showApertures && TargetPhotometryPanel.this.currentRadialProfile != null) {
                        n3 = (int)((double)TargetPhotometryPanel.this.currentRadialProfile.objectRadius * d);
                        graphics2D.setStroke(new BasicStroke(1.0f));
                        graphics2D.drawOval(n4 - n3, n5 - n3, n3 * 2, n3 * 2);
                        int n6 = (int)((double)TargetPhotometryPanel.this.currentRadialProfile.skyInnerRadius * d);
                        int n7 = (int)((double)TargetPhotometryPanel.this.currentRadialProfile.skyOuterRadius * d);
                        graphics2D.drawOval(n4 - n6, n5 - n6, n6 * 2, n6 * 2);
                        graphics2D.drawOval(n4 - n7, n5 - n7, n7 * 2, n7 * 2);
                    } else {
                        n3 = 4;
                        graphics2D.setStroke(new BasicStroke(2.0f));
                        graphics2D.drawOval(n4 - n3, n5 - n3, n3 * 2, n3 * 2);
                    }
                    graphics2D.setFont(new Font("SansSerif", 0, 10));
                    graphics2D.drawString(targetInfo.name, n4 + 6, n5 + 3);
                }
            }
            catch (Exception exception) {
                System.err.println("Error drawing variable stars: " + exception.getMessage());
            }
        }

        private void drawReferenceStars(Graphics2D graphics2D, int n, int n2, double d) {
            if (this.fitsHeader == null || this.stars == null) {
                return;
            }
            try {
                Double d2 = this.fitsHeader.getDoubleValue("CRVAL1");
                Double d3 = this.fitsHeader.getDoubleValue("CRVAL2");
                Double d4 = this.fitsHeader.getDoubleValue("CRPIX1");
                Double d5 = this.fitsHeader.getDoubleValue("CRPIX2");
                Double d6 = this.fitsHeader.getDoubleValue("CD1_1");
                Double d7 = this.fitsHeader.getDoubleValue("CD1_2");
                Double d8 = this.fitsHeader.getDoubleValue("CD2_1");
                Double d9 = this.fitsHeader.getDoubleValue("CD2_2");
                if (d2 == null || d3 == null || d4 == null || d5 == null) {
                    return;
                }
                if (d6 == null || d7 == null || d8 == null || d9 == null) {
                    Double d10 = this.fitsHeader.getDoubleValue("CDELT1");
                    Double object = this.fitsHeader.getDoubleValue("CDELT2");
                    Double d11 = this.fitsHeader.getDoubleValue("CROTA2");
                    if (d10 != null && object != null) {
                        double d12 = d11 != null ? Math.toRadians(d11) : 0.0;
                        d6 = d10 * Math.cos(d12);
                        d7 = -object.doubleValue() * Math.sin(d12);
                        d8 = d10 * Math.sin(d12);
                        d9 = object * Math.cos(d12);
                    } else {
                        return;
                    }
                }
                graphics2D.setColor(Color.RED);
                graphics2D.setStroke(new BasicStroke(2.0f));
                for (ReferenceStar referenceStar : this.stars) {
                    int n3;
                    double d13 = referenceStar.ra - d2;
                    double d14 = referenceStar.dec - d3;
                    d13 *= Math.cos(Math.toRadians(d3));
                    double d15 = d6 * d9 - d7 * d8;
                    if (Math.abs(d15) < 1.0E-10) continue;
                    double d16 = (d9 * d13 - d7 * d14) / d15;
                    double d17 = (-d8.doubleValue() * d13 + d6 * d14) / d15;
                    double d18 = d4 + d16 - 1.0;
                    double d19 = d5 + d17 - 1.0;
                    int n4 = n + (int)(d18 * d);
                    int n5 = n2 + (int)(d19 * d);
                    if (TargetPhotometryPanel.this.showApertures && TargetPhotometryPanel.this.currentRadialProfile != null) {
                        n3 = (int)((double)TargetPhotometryPanel.this.currentRadialProfile.objectRadius * d);
                        graphics2D.setStroke(new BasicStroke(1.0f));
                        graphics2D.drawOval(n4 - n3, n5 - n3, n3 * 2, n3 * 2);
                        int n6 = (int)((double)TargetPhotometryPanel.this.currentRadialProfile.skyInnerRadius * d);
                        int n7 = (int)((double)TargetPhotometryPanel.this.currentRadialProfile.skyOuterRadius * d);
                        graphics2D.drawOval(n4 - n6, n5 - n6, n6 * 2, n6 * 2);
                        graphics2D.drawOval(n4 - n7, n5 - n7, n7 * 2, n7 * 2);
                    } else {
                        n3 = 4;
                        graphics2D.setStroke(new BasicStroke(2.0f));
                        graphics2D.drawOval(n4 - n3, n5 - n3, n3 * 2, n3 * 2);
                    }
                    String string = referenceStar.label != null && !referenceStar.label.isEmpty() ? referenceStar.label : referenceStar.name;
                    graphics2D.setFont(new Font("SansSerif", 0, 10));
                    graphics2D.drawString(string, n4 + 6, n5 + 3);
                }
            }
            catch (Exception exception) {
                System.err.println("Error drawing reference stars: " + exception.getMessage());
            }
        }

        private void updateMouseInfo(int n, int n2) {
            if (TargetPhotometryPanel.this.mouseInfoLabel == null || TargetPhotometryPanel.this.currentImageData == null || this.image == null) {
                return;
            }
            int n3 = (int)((double)(n - this.imageOffsetX) / this.imageScale);
            int n4 = (int)((double)(n2 - this.imageOffsetY) / this.imageScale);
            if (n3 < 0 || n4 < 0 || n3 >= TargetPhotometryPanel.this.currentImageData[0].length || n4 >= TargetPhotometryPanel.this.currentImageData.length) {
                TargetPhotometryPanel.this.mouseInfoLabel.setText("x: -, y: -, intensity: -");
                return;
            }
            try {
                int n5 = TargetPhotometryPanel.this.currentImageData[n4][n3];
                TargetPhotometryPanel.this.mouseInfoLabel.setText(String.format("x: %d, y: %d, intensity: %d", n3, n4, n5));
            }
            catch (Exception exception) {
                TargetPhotometryPanel.this.mouseInfoLabel.setText("x: -, y: -, intensity: -");
            }
        }

        private void drawCoordinateBox(Graphics2D graphics2D, int n, int n2) {
            int n3 = (int)((double)(n - this.imageOffsetX) / this.imageScale);
            int n4 = (int)((double)(n2 - this.imageOffsetY) / this.imageScale);
            if (n3 < 0 || n4 < 0 || this.image == null || n3 >= this.image.getWidth() || n4 >= this.image.getHeight()) {
                return;
            }
            try {
                double d;
                Double d2 = this.fitsHeader.getDoubleValue("CRVAL1");
                Double d3 = this.fitsHeader.getDoubleValue("CRVAL2");
                Double d4 = this.fitsHeader.getDoubleValue("CRPIX1");
                Double d5 = this.fitsHeader.getDoubleValue("CRPIX2");
                Double d6 = this.fitsHeader.getDoubleValue("CD1_1");
                Double d7 = this.fitsHeader.getDoubleValue("CD1_2");
                Double d8 = this.fitsHeader.getDoubleValue("CD2_1");
                Double d9 = this.fitsHeader.getDoubleValue("CD2_2");
                if (d2 == null || d3 == null || d4 == null || d5 == null) {
                    return;
                }
                if (d6 == null || d7 == null || d8 == null || d9 == null) {
                    Double d10 = this.fitsHeader.getDoubleValue("CDELT1");
                    Double d11 = this.fitsHeader.getDoubleValue("CDELT2");
                    Double d12 = this.fitsHeader.getDoubleValue("CROTA2");
                    if (d10 != null && d11 != null) {
                        double d13 = d12 != null ? Math.toRadians(d12) : 0.0;
                        d6 = d10 * Math.cos(d13);
                        d7 = -d11.doubleValue() * Math.sin(d13);
                        d8 = d10 * Math.sin(d13);
                        d9 = d11 * Math.cos(d13);
                    } else {
                        return;
                    }
                }
                double d14 = (double)(n3 + 1) - d4;
                double d15 = (double)(n4 + 1) - d5;
                double d16 = d6 * d14 + d7 * d15;
                double d17 = d8 * d14 + d9 * d15;
                double d18 = d3 + d17;
                for (d = d2 + (d16 /= Math.cos(Math.toRadians(d3))); d < 0.0; d += 360.0) {
                }
                while (d >= 360.0) {
                    d -= 360.0;
                }
                String string = TargetPhotometryPanel.this.formatRA(d);
                String string2 = TargetPhotometryPanel.this.formatDec(d18);
                graphics2D.setFont(new Font("Monospaced", 0, 10));
                FontMetrics fontMetrics = graphics2D.getFontMetrics();
                String string3 = "RA: " + string + "  Dec: " + string2;
                int n5 = fontMetrics.stringWidth(string3);
                int n6 = fontMetrics.getHeight();
                int n7 = this.getWidth() - n5 - 15;
                int n8 = 10;
                int n9 = n5 + 10;
                int n10 = n6 + 6;
                graphics2D.setColor(new Color(0, 0, 0, 180));
                graphics2D.fillRect(n7, n8, n9, n10);
                graphics2D.setColor(Color.YELLOW);
                graphics2D.drawRect(n7, n8, n9, n10);
                graphics2D.drawString(string3, n7 + 5, n8 + n6);
            }
            catch (Exception exception) {
                // empty catch block
            }
        }
    }

    private static class TargetNightEntry {
        String target;
        LocalDate night;
        String status;
        int imageCount;

        TargetNightEntry(String string, LocalDate localDate, String string2, int n) {
            this.target = string;
            this.night = localDate;
            this.status = string2;
            this.imageCount = n;
        }
    }

    private static class PhotometryResult {
        String imageId;
        String name;
        String type;
        String filter;
        double catalogRA;
        double catalogDec;
        double centroidX;
        double centroidY;
        double centroidRA;
        double centroidDec;
        double starFlux;
        double skyFlux;
        double skyPerPixel;
        double netFlux;
        double instMag;
        double magError;
        double calibratedMag;
        double zeropoint;
        double zpUncertainty;

        PhotometryResult(String string, String string2, String string3, String string4, double d, double d2, double d3, double d4, double d5, double d6, double d7, double d8, double d9, double d10, double d11, double d12) {
            this.imageId = string;
            this.name = string2;
            this.type = string3;
            this.filter = string4;
            this.catalogRA = d;
            this.catalogDec = d2;
            this.centroidX = d3;
            this.centroidY = d4;
            this.centroidRA = d5;
            this.centroidDec = d6;
            this.starFlux = d7;
            this.skyFlux = d8;
            this.skyPerPixel = d9;
            this.netFlux = d10;
            this.instMag = d11;
            this.magError = d12;
            this.calibratedMag = Double.NaN;
            this.zeropoint = Double.NaN;
            this.zpUncertainty = Double.NaN;
        }

        void applyZeropoint(double d, double d2) {
            this.zeropoint = d;
            this.zpUncertainty = d2;
            this.calibratedMag = this.instMag + d;
        }
    }

    private static class LightcurvePoint {
        double hjd;
        double mag;
        double error;

        LightcurvePoint(double d, double d2, double d3) {
            this.hjd = d;
            this.mag = d2;
            this.error = d3;
        }
    }

    private static class LightcurvePlotPanel
    extends JPanel {
        private String varName;
        private List<LightcurvePoint> varPoints;
        private List<LightcurvePoint> checkPoints;

        LightcurvePlotPanel(String string, List<LightcurvePoint> list, List<LightcurvePoint> list2) {
            this.varName = string;
            this.varPoints = list;
            this.checkPoints = list2;
            this.setBackground(Color.WHITE);
        }

        @Override
        protected void paintComponent(Graphics graphics) {
            int n;
            int n2;
            int n3;
            int n4;
            double d;
            int n5;
            super.paintComponent(graphics);
            Graphics2D graphics2D = (Graphics2D)graphics;
            graphics2D.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
            if (this.varPoints.isEmpty()) {
                return;
            }
            int n6 = this.getWidth();
            int n7 = this.getHeight();
            int n8 = 60;
            int n9 = n6 - 2 * n8;
            int n10 = n7 - 2 * n8;
            double d2 = Double.MAX_VALUE;
            double d3 = Double.MIN_VALUE;
            double d4 = Double.MAX_VALUE;
            double d5 = Double.MIN_VALUE;
            for (LightcurvePoint lightcurvePoint : this.varPoints) {
                d2 = Math.min(d2, lightcurvePoint.hjd);
                d3 = Math.max(d3, lightcurvePoint.hjd);
                d4 = Math.min(d4, lightcurvePoint.mag - lightcurvePoint.error);
                d5 = Math.max(d5, lightcurvePoint.mag + lightcurvePoint.error);
            }
            for (LightcurvePoint lightcurvePoint : this.checkPoints) {
                d2 = Math.min(d2, lightcurvePoint.hjd);
                d3 = Math.max(d3, lightcurvePoint.hjd);
                d4 = Math.min(d4, lightcurvePoint.mag - lightcurvePoint.error);
                d5 = Math.max(d5, lightcurvePoint.mag + lightcurvePoint.error);
            }
            double d6 = d3 - d2;
            if (d6 < 0.1) {
                d6 = 0.1;
            }
            d2 -= d6 * 0.05;
            d3 += d6 * 0.05;
            double d7 = d5 - d4;
            if (d7 < 0.1) {
                d7 = 0.1;
            }
            d4 -= d7 * 0.1;
            d5 += d7 * 0.1;
            graphics2D.setColor(Color.BLACK);
            graphics2D.setStroke(new BasicStroke(2.0f));
            graphics2D.drawLine(n8, n8, n8, n7 - n8);
            graphics2D.drawLine(n8, n7 - n8, n6 - n8, n7 - n8);
            graphics2D.setFont(new Font("SansSerif", 1, 14));
            FontMetrics fontMetrics = graphics2D.getFontMetrics();
            String string = "HJD - " + String.format("%.0f", Math.floor(d2));
            int n11 = fontMetrics.stringWidth(string);
            graphics2D.drawString(string, (n6 - n11) / 2, n7 - 10);
            String string2 = "Magnitude";
            graphics2D.rotate(-1.5707963267948966);
            graphics2D.drawString(string2, -(n7 + fontMetrics.stringWidth(string2)) / 2, 20);
            graphics2D.rotate(1.5707963267948966);
            graphics2D.setFont(new Font("SansSerif", 1, 16));
            String string3 = "Lightcurve: " + this.varName;
            int n12 = graphics2D.getFontMetrics().stringWidth(string3);
            graphics2D.drawString(string3, (n6 - n12) / 2, 25);
            graphics2D.setColor(Color.LIGHT_GRAY);
            graphics2D.setStroke(new BasicStroke(1.0f));
            graphics2D.setFont(new Font("SansSerif", 0, 10));
            int n13 = 5;
            for (n5 = 0; n5 <= n13; ++n5) {
                d = d4 + (d5 - d4) * (double)n5 / (double)n13;
                n4 = n8 + n10 - (int)((d - d4) / (d5 - d4) * (double)n10);
                graphics2D.setColor(Color.LIGHT_GRAY);
                graphics2D.drawLine(n8, n4, n6 - n8, n4);
                graphics2D.setColor(Color.BLACK);
                graphics2D.drawLine(n8 - 5, n4, n8, n4);
                String string4 = String.format("%.2f", d);
                graphics2D.drawString(string4, n8 - 45, n4 + 4);
            }
            n5 = 5;
            d = Math.floor(d2);
            for (n4 = 0; n4 <= n5; ++n4) {
                double d8 = d2 + (d3 - d2) * (double)n4 / (double)n5;
                n3 = n8 + (int)((d8 - d2) / (d3 - d2) * (double)n9);
                graphics2D.setColor(Color.LIGHT_GRAY);
                graphics2D.drawLine(n3, n8, n3, n7 - n8);
                graphics2D.setColor(Color.BLACK);
                graphics2D.drawLine(n3, n7 - n8, n3, n7 - n8 + 5);
                String string5 = String.format("%.3f", d8 - d);
                n2 = graphics2D.getFontMetrics().stringWidth(string5);
                graphics2D.drawString(string5, n3 - n2 / 2, n7 - n8 + 20);
            }
            if (!this.checkPoints.isEmpty()) {
                graphics2D.setColor(Color.GRAY);
                for (LightcurvePoint lightcurvePoint : this.checkPoints) {
                    n = n8 + (int)((lightcurvePoint.hjd - d2) / (d3 - d2) * (double)n9);
                    n3 = n8 + n10 - (int)((lightcurvePoint.mag - d4) / (d5 - d4) * (double)n10);
                    int n14 = n8 + n10 - (int)((lightcurvePoint.mag - lightcurvePoint.error - d4) / (d5 - d4) * (double)n10);
                    n2 = n8 + n10 - (int)((lightcurvePoint.mag + lightcurvePoint.error - d4) / (d5 - d4) * (double)n10);
                    graphics2D.drawLine(n, n14, n, n2);
                    graphics2D.drawLine(n - 3, n14, n + 3, n14);
                    graphics2D.drawLine(n - 3, n2, n + 3, n2);
                    graphics2D.fillOval(n - 3, n3 - 3, 6, 6);
                }
            }
            graphics2D.setColor(Color.BLUE);
            for (LightcurvePoint lightcurvePoint : this.varPoints) {
                n = n8 + (int)((lightcurvePoint.hjd - d2) / (d3 - d2) * (double)n9);
                n3 = n8 + n10 - (int)((lightcurvePoint.mag - d4) / (d5 - d4) * (double)n10);
                int n15 = n8 + n10 - (int)((lightcurvePoint.mag - lightcurvePoint.error - d4) / (d5 - d4) * (double)n10);
                n2 = n8 + n10 - (int)((lightcurvePoint.mag + lightcurvePoint.error - d4) / (d5 - d4) * (double)n10);
                graphics2D.drawLine(n, n15, n, n2);
                graphics2D.drawLine(n - 3, n15, n + 3, n15);
                graphics2D.drawLine(n - 3, n2, n + 3, n2);
                graphics2D.fillOval(n - 4, n3 - 4, 8, 8);
            }
            int n16 = n6 - n8 - 120;
            int n17 = n8 + 20;
            graphics2D.setColor(Color.BLUE);
            graphics2D.fillOval(n16, n17 - 4, 8, 8);
            graphics2D.setColor(Color.BLACK);
            graphics2D.setFont(new Font("SansSerif", 0, 12));
            graphics2D.drawString(this.varName, n16 + 15, n17 + 4);
            if (!this.checkPoints.isEmpty()) {
                graphics2D.setColor(Color.GRAY);
                graphics2D.fillOval(n16, n17 + 16, 6, 6);
                graphics2D.setColor(Color.BLACK);
                graphics2D.drawString("CHECK star", n16 + 15, n17 + 24);
            }
        }
    }

    private static class RadialProfilePlot
    extends JPanel {
        private RadialProfile profile;

        RadialProfilePlot(RadialProfile radialProfile) {
            this.profile = radialProfile;
            this.setBackground(Color.WHITE);
        }

        void updateProfile(RadialProfile radialProfile) {
            this.profile = radialProfile;
            this.repaint();
        }

        @Override
        protected void paintComponent(Graphics graphics) {
            int n;
            int n2;
            int n3;
            super.paintComponent(graphics);
            Graphics2D graphics2D = (Graphics2D)graphics;
            graphics2D.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
            graphics2D.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_ON);
            int n4 = this.getWidth();
            int n5 = this.getHeight();
            int n6 = 80;
            int n7 = 40;
            int n8 = 40;
            int n9 = 60;
            int n10 = n4 - n6 - n7;
            int n11 = n5 - n8 - n9;
            double d = 0.0;
            double[] dArray = this.profile.intensity;
            int n12 = dArray.length;
            for (n3 = 0; n3 < n12; ++n3) {
                double d2 = dArray[n3];
                d = Math.max(d, d2);
            }
            double d3 = Math.pow(10.0, Math.floor(Math.log10(d)));
            d = Math.ceil(d / d3) * d3;
            graphics2D.setColor(Color.WHITE);
            graphics2D.fillRect(n6, n8, n10, n11);
            graphics2D.setColor(new Color(230, 230, 230));
            graphics2D.setStroke(new BasicStroke(1.0f));
            for (n3 = 0; n3 <= 10; ++n3) {
                int n13 = n5 - n9 - n3 * n11 / 10;
                graphics2D.drawLine(n6, n13, n6 + n10, n13);
                int n14 = n6 + n3 * n10 / 10;
                graphics2D.drawLine(n14, n8, n14, n5 - n9);
            }
            graphics2D.setColor(Color.BLACK);
            graphics2D.setStroke(new BasicStroke(2.0f));
            graphics2D.drawLine(n6, n5 - n9, n6 + n10, n5 - n9);
            graphics2D.drawLine(n6, n5 - n9, n6, n8);
            graphics2D.setFont(new Font("SansSerif", 1, 13));
            graphics2D.drawString("Radius (pixels)", n4 / 2 - 50, n5 - 10);
            graphics2D.rotate(-1.5707963267948966);
            graphics2D.drawString("Mean Intensity (ADU)", -n5 / 2 - 70, 25);
            graphics2D.rotate(1.5707963267948966);
            graphics2D.setFont(new Font("SansSerif", 0, 11));
            for (n3 = 0; n3 <= 10; ++n3) {
                int n15 = n5 - n9 - n3 * n11 / 10;
                graphics2D.setColor(Color.BLACK);
                graphics2D.setStroke(new BasicStroke(2.0f));
                graphics2D.drawLine(n6 - 5, n15, n6, n15);
                String string = String.format("%.0f", d * (double)n3 / 10.0);
                int n16 = graphics2D.getFontMetrics().stringWidth(string);
                graphics2D.drawString(string, n6 - n16 - 10, n15 + 4);
            }
            n3 = this.profile.intensity.length - 1;
            for (n2 = 0; n2 <= 10; ++n2) {
                int n17 = n6 + n2 * n10 / 10;
                graphics2D.setColor(Color.BLACK);
                graphics2D.setStroke(new BasicStroke(2.0f));
                graphics2D.drawLine(n17, n5 - n9, n17, n5 - n9 + 5);
                String string = String.format("%d", n3 * n2 / 10);
                n = graphics2D.getFontMetrics().stringWidth(string);
                graphics2D.drawString(string, n17 - n / 2, n5 - n9 + 20);
            }
            graphics2D.setColor(new Color(0, 102, 204));
            graphics2D.setStroke(new BasicStroke(1.5f));
            for (n2 = 0; n2 < this.profile.intensity.length - 1; ++n2) {
                int n18 = n6 + n2 * n10 / n3;
                int n19 = n5 - n9 - (int)(this.profile.intensity[n2] * (double)n11 / d);
                n = n6 + (n2 + 1) * n10 / n3;
                int n20 = n5 - n9 - (int)(this.profile.intensity[n2 + 1] * (double)n11 / d);
                graphics2D.drawLine(n18, n19, n, n20);
            }
            graphics2D.setColor(new Color(0, 102, 204));
            for (n2 = 0; n2 < this.profile.intensity.length; ++n2) {
                int n21 = n6 + n2 * n10 / n3;
                int n22 = n5 - n9 - (int)(this.profile.intensity[n2] * (double)n11 / d);
                graphics2D.fillOval(n21 - 3, n22 - 3, 6, 6);
            }
            graphics2D.setColor(new Color(0, 150, 0));
            graphics2D.setStroke(new BasicStroke(2.0f, 0, 0, 10.0f, new float[]{8.0f, 4.0f}, 0.0f));
            n2 = n6 + this.profile.objectRadius * n10 / n3;
            if (n2 > n6 && n2 < n6 + n10) {
                graphics2D.drawLine(n2, n5 - n9, n2, n8);
                graphics2D.setFont(new Font("SansSerif", 1, 11));
                graphics2D.drawString("Obj", n2 + 5, n8 + 15);
            }
            graphics2D.setColor(new Color(204, 102, 0));
            int n23 = n6 + this.profile.skyInnerRadius * n10 / n3;
            int n24 = n6 + this.profile.skyOuterRadius * n10 / n3;
            if (n23 > n6 && n23 < n6 + n10) {
                graphics2D.drawLine(n23, n5 - n9, n23, n8);
            }
            if (n24 > n6 && n24 < n6 + n10) {
                graphics2D.drawLine(n24, n5 - n9, n24, n8);
            }
            if (n23 > n6 && n24 < n6 + n10) {
                n = (n23 + n24) / 2;
                graphics2D.setFont(new Font("SansSerif", 1, 11));
                graphics2D.drawString("Sky", n - 12, n8 + 15);
            }
            graphics2D.setColor(Color.BLACK);
            graphics2D.setStroke(new BasicStroke(2.0f));
            graphics2D.drawRect(n6, n8, n10, n11);
        }
    }

    private static class TargetInfo {
        String name;
        String auid;
        double ra;
        double dec;
        String variabilityType;
        double period;
        double epoch;
        String maxMag;
        String minMag;
        String constellation;

        TargetInfo(String string) {
            this.name = string;
        }
    }

    private static class FluxMeasurement {
        double starFlux;
        double skyFlux;
        double skyPerPixel;
        double netFlux;
        double error;

        FluxMeasurement(double d, double d2, double d3, double d4, double d5) {
            this.starFlux = d;
            this.skyFlux = d2;
            this.skyPerPixel = d3;
            this.netFlux = d4;
            this.error = d5;
        }
    }
}

