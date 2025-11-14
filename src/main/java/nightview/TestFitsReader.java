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

import java.io.File;
import nom.tam.fits.Fits;
import nom.tam.fits.Header;
import nom.tam.fits.ImageData;
import nom.tam.fits.ImageHDU;

public class TestFitsReader {
    public static void main(String[] stringArray) {
        if (stringArray.length < 1) {
            System.err.println("Usage: java TestFitsReader <fits_file>");
            System.exit(1);
        }
        try {
            File file = new File(stringArray[0]);
            System.out.println("Reading: " + file.getName());
            Fits fits = new Fits(file);
            ImageHDU imageHDU = (ImageHDU)fits.getHDU(0);
            Object object = ((ImageData)imageHDU.getData()).getData();
            Header header = imageHDU.getHeader();
            System.out.println("BITPIX: " + header.getIntValue("BITPIX"));
            System.out.println("NAXIS1: " + header.getIntValue("NAXIS1"));
            System.out.println("NAXIS2: " + header.getIntValue("NAXIS2"));
            System.out.println("Data class: " + object.getClass().getName());
            int[][] nArray = TestFitsReader.convertToIntArray(object);
            int n = Integer.MAX_VALUE;
            int n2 = Integer.MIN_VALUE;
            int n3 = 0;
            int n4 = 0;
            for (int i = 0; i < nArray.length; ++i) {
                for (int j = 0; j < nArray[i].length; ++j) {
                    if (nArray[i][j] < n) {
                        n = nArray[i][j];
                    }
                    if (nArray[i][j] <= n2) continue;
                    n2 = nArray[i][j];
                    n3 = i;
                    n4 = j;
                }
            }
            System.out.println("\nConverted to int[][]:");
            System.out.println("Min: " + n);
            System.out.println("Max: " + n2);
            System.out.println("Max at (y=" + n3 + ", x=" + n4 + ")");
            System.out.println("\nSample pixel values:");
            System.out.println("At (2428, 2402): " + nArray[2428][2402]);
            System.out.println("At (455, 1201): " + nArray[455][1201]);
            System.out.println("At (2166, 3135): " + nArray[2166][3135]);
            fits.close();
        }
        catch (Exception exception) {
            exception.printStackTrace();
        }
    }

    private static int[][] convertToIntArray(Object object) {
        if (object instanceof short[][]) {
            short[][] sArray = (short[][])object;
            int[][] nArray = new int[sArray.length][sArray[0].length];
            for (int i = 0; i < sArray.length; ++i) {
                for (int j = 0; j < sArray[i].length; ++j) {
                    nArray[i][j] = sArray[i][j] & 0xFFFF;
                }
            }
            System.out.println("Converted from short[][]");
            return nArray;
        }
        if (object instanceof int[][]) {
            System.out.println("Already int[][]");
            return (int[][])object;
        }
        if (object instanceof float[][]) {
            float[][] fArray = (float[][])object;
            int[][] nArray = new int[fArray.length][fArray[0].length];
            float f = Float.MAX_VALUE;
            float f2 = Float.MIN_VALUE;
            for (int i = 0; i < fArray.length; ++i) {
                for (int j = 0; j < fArray[i].length; ++j) {
                    if (fArray[i][j] < f) {
                        f = fArray[i][j];
                    }
                    if (fArray[i][j] > f2) {
                        f2 = fArray[i][j];
                    }
                    nArray[i][j] = (int)fArray[i][j];
                }
            }
            System.out.println("Converted from float[][] - float min=" + f + ", max=" + f2);
            return nArray;
        }
        throw new IllegalArgumentException("Unsupported data type: " + String.valueOf(object.getClass()));
    }
}

