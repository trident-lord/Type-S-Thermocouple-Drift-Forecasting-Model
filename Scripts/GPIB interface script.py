import pyvisa
import time
import csv
from datetime import datetime

def log_fluke_8508a_emf_continuous(output_csv='fluke8508a_log.csv', gpib_address='GPIB0::22::INSTR',
                                    interval_sec=1, save_every_sec=60):
    import os

    rm = pyvisa.ResourceManager()
    meter = rm.open_resource(gpib_address)

    # Setup Fluke 8508A for DC voltage in high resolution
    meter.write("FUNC 'VOLT:DC'")
    meter.write("VOLT:DC:RANG 0.1")
    meter.write("VOLT:DC:RES 0.00000001")

    buffer = []
    last_save_time = time.time()

    # If file doesn't exist, write header
    if not os.path.exists(output_csv):
        with open(output_csv, mode='w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(["Timestamp", "EMF (mV)"])

    print("[⚡] Logging started. Press Ctrl+C to stop.")

    try:
        while True:
            try:
                emf_v = float(meter.query("READ?").strip())
                emf_mv = emf_v * 1000
                timestamp = datetime.now().isoformat()
                buffer.append([timestamp, emf_mv])
                print(f"{timestamp} - EMF: {emf_mv:.8f} mV")
            except Exception as e:
                print(f"[⚠] Read error: {e}")

            time.sleep(interval_sec)

            # Autosave every N seconds
            if time.time() - last_save_time >= save_every_sec and buffer:
                with open(output_csv, mode='a', newline='') as f:
                    writer = csv.writer(f)
                    writer.writerows(buffer)
                print(f"[💾] Autosaved {len(buffer)} rows to {output_csv}")
                buffer.clear()
                last_save_time = time.time()

    except KeyboardInterrupt:
        print("\n[🛑] Logging stopped by user.")
        if buffer:
            with open(output_csv, mode='a', newline='') as f:
                writer = csv.writer(f)
                writer.writerows(buffer)
            print(f"[💾] Final autosave complete: {len(buffer)} rows.")
        print(f"[✅] Log file saved to: {output_csv}")


log_fluke_8508a_emf_continuous(output_csv='8508A_al_point_log.csv')
