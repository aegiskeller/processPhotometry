#!/usr/bin/env python3
"""Quick diagnostic to check pixel values in FITS files."""

import sys
import struct

def check_fits_file(filename):
    """Check min/max pixel values in a FITS file without external dependencies."""
    print(f"\nChecking: {filename}")
    
    with open(filename, 'rb') as f:
        # Read FITS header (2880 bytes per block)
        header_blocks = []
        while True:
            block = f.read(2880)
            if not block:
                break
            header_blocks.append(block)
            if b'END' in block:
                # Find END card and stop reading header
                end_pos = block.find(b'END')
                if end_pos >= 0 and (end_pos % 80 == 0):
                    break
        
        # Parse key header keywords
        header = b''.join(header_blocks)
        
        def get_keyword(kw):
            pattern = kw.encode() + b' '
            idx = header.find(pattern)
            if idx < 0:
                return None
            card = header[idx:idx+80].decode('ascii', errors='ignore')
            if '=' not in card:
                return None
            value = card.split('=')[1].split('/')[0].strip()
            try:
                return int(value) if value.isdigit() or (value[0] == '-' and value[1:].isdigit()) else value.strip("'\" ")
            except:
                return value.strip("'\" ")
        
        bitpix = get_keyword('BITPIX')
        naxis1 = get_keyword('NAXIS1')
        naxis2 = get_keyword('NAXIS2')
        
        print(f"  BITPIX: {bitpix}")
        print(f"  NAXIS1: {naxis1}, NAXIS2: {naxis2}")
        
        if bitpix and naxis1 and naxis2:
            bytes_per_pixel = abs(bitpix) // 8
            total_pixels = naxis1 * naxis2
            
            # Read a sample of pixel data
            data_start = len(b''.join(header_blocks))
            f.seek(data_start)
            
            # Read first 1000 pixels and last 1000 pixels as samples
            sample_size = min(1000, total_pixels)
            
            if bitpix == 16:  # 16-bit signed
                fmt = '>h'  # big-endian short
                bytes_per_pixel = 2
            elif bitpix == -32:  # 32-bit float
                fmt = '>f'
                bytes_per_pixel = 4
            else:
                print(f"  Unsupported BITPIX: {bitpix}")
                return
            
            # Read samples (read more samples distributed across image)
            samples = []
            step = max(1, total_pixels // 10000)  # Sample ~10000 pixels across image
            for i in range(0, min(total_pixels, 100000), step):
                f.seek(data_start + i * bytes_per_pixel)
                data = f.read(bytes_per_pixel)
                if not data or len(data) < bytes_per_pixel:
                    break
                value = struct.unpack(fmt, data)[0]
                if bitpix == 16:
                    # Convert to unsigned if needed
                    value = value if value >= 0 else 65536 + value
                samples.append(value)
            
            if samples:
                print(f"  Sample of {len(samples)} pixels:")
                print(f"    Min: {min(samples)}")
                print(f"    Max: {max(samples)}")
                print(f"    Mean: {sum(samples)/len(samples):.1f}")

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: check_fits_values.py <fits_file> [<fits_file> ...]")
        sys.exit(1)
    
    for filename in sys.argv[1:]:
        try:
            check_fits_file(filename)
        except Exception as e:
            print(f"Error reading {filename}: {e}")
