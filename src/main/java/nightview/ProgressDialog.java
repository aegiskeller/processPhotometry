/*
 * Decompiled with CFR 0.152.
 */
package nightview;

import java.awt.BorderLayout;
import java.awt.Component;
import java.awt.Font;
import java.awt.Frame;
import javax.swing.BorderFactory;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import javax.swing.SwingUtilities;

public class ProgressDialog
extends JDialog {
    private JLabel statusLabel;
    private JTextArea logArea;
    private JScrollPane scrollPane;

    public ProgressDialog(Frame frame, String string) {
        super(frame, string, false);
        this.initializeUI();
    }

    private void initializeUI() {
        this.setLayout(new BorderLayout());
        this.statusLabel = new JLabel("Initializing...", 0);
        this.statusLabel.setBorder(BorderFactory.createEmptyBorder(10, 10, 10, 10));
        this.statusLabel.setFont(this.statusLabel.getFont().deriveFont(1));
        this.add((Component)this.statusLabel, "North");
        this.logArea = new JTextArea(15, 50);
        this.logArea.setEditable(false);
        this.logArea.setFont(new Font("Monospaced", 0, 12));
        this.scrollPane = new JScrollPane(this.logArea);
        this.scrollPane.setVerticalScrollBarPolicy(22);
        this.add((Component)this.scrollPane, "Center");
        this.setDefaultCloseOperation(0);
        this.pack();
        this.setLocationRelativeTo(this.getParent());
        this.setResizable(true);
    }

    public void updateProgress(String string) {
        SwingUtilities.invokeLater(() -> {
            this.statusLabel.setText(string);
            this.logArea.append(string + "\n");
            this.logArea.setCaretPosition(this.logArea.getDocument().getLength());
        });
    }

    public void clearLog() {
        SwingUtilities.invokeLater(() -> this.logArea.setText(""));
    }
}

