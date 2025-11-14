/*
 * Decompiled with CFR 0.152.
 */
package nightview;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.GridLayout;
import java.time.LocalDate;
import java.time.YearMonth;
import java.time.format.DateTimeFormatter;
import java.util.HashSet;
import java.util.Set;
import javax.swing.JButton;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.border.LineBorder;
import nightview.NightViewApp;

public class NightCalendarPanel
extends JPanel {
    private NightViewApp parentApp;
    private Set<LocalDate> availableNights;
    private Set<LocalDate> processedNights;
    private LocalDate selectedDate;
    private LocalDate currentMonth;
    private JLabel monthLabel;
    private JPanel calendarGrid;

    public NightCalendarPanel(NightViewApp nightViewApp) {
        this.parentApp = nightViewApp;
        this.availableNights = new HashSet<LocalDate>();
        this.processedNights = new HashSet<LocalDate>();
        this.currentMonth = LocalDate.now().withDayOfMonth(1);
        this.initializeUI();
    }

    private void initializeUI() {
        this.setLayout(new BorderLayout());
        JPanel jPanel = this.createNavigationPanel();
        this.add((Component)jPanel, "North");
        this.calendarGrid = new JPanel(new GridLayout(7, 7, 1, 1));
        this.calendarGrid.setBackground(Color.LIGHT_GRAY);
        this.add((Component)this.calendarGrid, "Center");
        this.updateCalendarDisplay();
    }

    private JPanel createNavigationPanel() {
        JPanel jPanel = new JPanel(new BorderLayout());
        JButton jButton = new JButton("<");
        jButton.addActionListener(actionEvent -> {
            this.currentMonth = this.currentMonth.minusMonths(1L);
            this.updateCalendarDisplay();
        });
        JButton jButton2 = new JButton(">");
        jButton2.addActionListener(actionEvent -> {
            this.currentMonth = this.currentMonth.plusMonths(1L);
            this.updateCalendarDisplay();
        });
        this.monthLabel = new JLabel("", 0);
        this.monthLabel.setFont(this.monthLabel.getFont().deriveFont(1, 14.0f));
        jPanel.add((Component)jButton, "West");
        jPanel.add((Component)this.monthLabel, "Center");
        jPanel.add((Component)jButton2, "East");
        return jPanel;
    }

    private void updateCalendarDisplay() {
        int n;
        int n2;
        this.calendarGrid.removeAll();
        this.monthLabel.setText(this.currentMonth.format(DateTimeFormatter.ofPattern("MMMM yyyy")));
        String[] stringArray = new String[]{"Sun", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat"};
        for (String string : stringArray) {
            JLabel jLabel = new JLabel(string, 0);
            jLabel.setFont(jLabel.getFont().deriveFont(1));
            jLabel.setBorder(new LineBorder(Color.GRAY));
            jLabel.setOpaque(true);
            jLabel.setBackground(Color.LIGHT_GRAY);
            this.calendarGrid.add(jLabel);
        }
        YearMonth yearMonth = YearMonth.from(this.currentMonth);
        LocalDate localDate = yearMonth.atDay(1);
        int n3 = localDate.getDayOfWeek().getValue() % 7;
        for (n2 = 0; n2 < n3; ++n2) {
            this.calendarGrid.add(this.createEmptyDayButton());
        }
        n2 = yearMonth.lengthOfMonth();
        for (n = 1; n <= n2; ++n) {
            LocalDate localDate2 = localDate.withDayOfMonth(n);
            this.calendarGrid.add(this.createDayButton(localDate2));
        }
        for (n = this.calendarGrid.getComponentCount(); n < 49; ++n) {
            this.calendarGrid.add(this.createEmptyDayButton());
        }
        this.calendarGrid.revalidate();
        this.calendarGrid.repaint();
    }

    private JButton createDayButton(LocalDate localDate) {
        JButton jButton = new JButton(String.valueOf(localDate.getDayOfMonth()));
        jButton.setPreferredSize(new Dimension(40, 30));
        if (this.processedNights.contains(localDate)) {
            jButton.setBackground(new Color(70, 130, 180));
            jButton.setForeground(Color.WHITE);
            jButton.setOpaque(true);
            jButton.setBorder(new LineBorder(new Color(25, 25, 112), 2));
            jButton.setFont(jButton.getFont().deriveFont(1));
            jButton.setToolTipText("FULLY PROCESSED - " + localDate.format(DateTimeFormatter.ofPattern("yyyy-MM-dd")));
        } else if (this.availableNights.contains(localDate)) {
            jButton.setBackground(Color.GREEN);
            jButton.setOpaque(true);
            jButton.setBorder(new LineBorder(Color.GREEN.darker()));
            jButton.setToolTipText("Night has data available");
        } else {
            jButton.setBackground(null);
            jButton.setOpaque(false);
            jButton.setBorder(new LineBorder(Color.GRAY));
            jButton.setToolTipText("No data available - click to select anyway");
        }
        if (localDate.equals(this.selectedDate)) {
            jButton.setBorder(new LineBorder(Color.BLUE, 2));
        }
        jButton.addActionListener(actionEvent -> {
            this.selectedDate = localDate;
            this.parentApp.onNightSelected(localDate);
            this.updateCalendarDisplay();
        });
        jButton.setEnabled(true);
        return jButton;
    }

    private JButton createEmptyDayButton() {
        JButton jButton = new JButton();
        jButton.setEnabled(false);
        jButton.setVisible(false);
        jButton.setPreferredSize(new Dimension(40, 30));
        return jButton;
    }

    public void setAvailableNights(Set<LocalDate> set) {
        this.availableNights = set;
        this.updateCalendarDisplay();
    }

    public void setProcessedNights(Set<LocalDate> set) {
        this.processedNights = set;
        this.updateCalendarDisplay();
    }

    public LocalDate getSelectedDate() {
        return this.selectedDate;
    }
}

