import sqlite3

conn = sqlite3.connect('wombatpipeline.db')
cursor = conn.cursor()

# Check if chart_id column exists
cursor.execute("PRAGMA table_info(target)")
columns = [row[1] for row in cursor.fetchall()]

if 'chart_id' not in columns:
    print("Adding chart_id column to target table...")
    cursor.execute("ALTER TABLE target ADD COLUMN chart_id TEXT")
    conn.commit()
    print("✓ Column added")
else:
    print("chart_id column already exists")

# Reset apphot status for RW PsA
cursor.execute("UPDATE target SET apphot = 0 WHERE date='2025-09-30' AND target_name='RW PsA'")
conn.commit()
print("✓ Reset RW PsA status")

# Check result
cursor.execute("SELECT date, target_name, filter, apphot, chart_id FROM target WHERE date='2025-09-30'")
row = cursor.fetchone()
if row:
    print(f"\n{row[0]} | {row[1]} | {row[2]} | apphot={row[3]} | chart_id={row[4]}")

conn.close()
