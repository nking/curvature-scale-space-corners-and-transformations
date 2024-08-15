package algorithms.imageProcessing.transform;

import algorithms.imageProcessing.ImageExt;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.matrix.MatrixUtil;
import algorithms.util.FormatArray;
import algorithms.util.ResourceFinder;
import no.uib.cipr.matrix.NotConvergedException;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;

/**
 * convenience methods for use with some of MASKS book data from examples_code.
 * Ma, Soatto, Kosecka,& Sastry (MASKS)
 * "An Invitation to Computer Vision, From Images to Geometric Models"
 *
 *the oldhouse images have rotation and translation
 *the yos images have zoom (translation in z direction, perpendicular to image  plane)
 * @author nichole
 */
public class MASKSData {
    
    public static String MASKSDIR = "masks";
    protected final static String sep = System.getProperty("file.separator");

    public static double[][] getFMOldhouse2A2_0_84() {
        double[][] fm = new double[3][3];
        fm[0] = new double[]{1.83300372e-06, -3.50706415e-04,  1.12755815e-01};
        fm[1] = new double[]{3.33245228e-04,  1.47091793e-05,  3.60295158e-01};
        fm[2] = new double[]{-1.01657683e-01, -3.55225108e-01,  8.49090177e-01};
        return fm;
    }

    public static double[][] getXOldhouse2A2000000() throws IOException {
        double[][] x = new double[3][72];
        x[0] = new double[]{537, 535, 636, 561, 631, 569, 502, 415, 343, 408, 518, 428,  99,
                610, 406, 504, 610, 400, 533, 491, 504, 232, 454, 609, 219, 441,
                617,  69, 119, 238, 276, 275, 563, 139, 108, 671, 188, 608, 341,
                605,  65, 634, 282, 228, 258,  70, 453, 254, 335, 198, 465, 213,
                496, 223, 165, 467, 181, 175, 311, 488, 231, 174, 150, 347, 642,
                220, 140, 267, 218, 182, 222, 670};
        x[1] = new double[]{283, 325, 340,  72, 267, 267, 255, 291, 319, 326, 239, 207, 429,
                285, 258, 148, 256, 229, 394, 302, 172,  71, 298, 312, 185, 241,
                412, 399, 418, 147, 420, 162, 206, 408, 393, 300,  52, 341, 271,
                376, 421, 314, 187,  94,  85,  50, 210, 249, 354, 168, 332,  80,
                373, 208, 165, 192, 403, 123, 203, 208, 429,  88, 244, 247, 169,
                383, 200, 120, 237, 289, 117, 188};
        x[2] = new double[x[0].length];
        Arrays.fill(x[2], 1);

        return x;
    }
    public static double[][] getXOldhouse2A2000084() throws IOException {
        double[][] x = new double[3][72];
        x[0] = new double[]{509.80793275, 508.22116716, 605.40931421, 532.74093715,
                599.70591021, 539.85224055, 475.94369746, 393.71456081,
                323.6288516 , 387.45697441, 490.91459628, 405.27798391,
                77.04621854, 579.1281747 , 384.89933471, 477.18767731,
                578.992311  , 378.87534753, 506.37131362, 465.99034968,
                477.35854666, 204.57073232, 430.75375104, 578.25400554,
                192.60865461, 417.80559304, 586.0018212 ,  50.03976583,
                96.02938843, 211.57280788, 247.32879836, 250.04704712,
                533.76721182, 115.32016067,  86.56217703, 637.99638315,
                162.10222369, 577.33520137, 321.01842593, 574.55189962,
                45.38463559, 602.81163809, 257.65541669, 200.6533172 ,
                231.56594661,  50.21673813, 428.97906425, 229.03610132,
                315.80943251, 172.63069221, 441.4786747 , 185.43935249,
                471.24502569, 196.93614127, 142.35341896, 442.15778526,
                157.20138133, 150.83363672, 288.13402747, 462.20995809,
                205.31541276, 149.53911958, 129.75021279, 326.58452264,
                608.87586674, 196.32518647, 119.90156891, 241.27646011,
                192.27041858, 159.36103496, 194.71790821, 635.77668323};
        x[1] = new double[]{269.86086148, 311.08973843, 325.00402395,  62.95771076,
                253.95860842, 254.30166208, 243.05461433, 278.73276534,
                306.30010944, 312.80337504, 227.30171578, 196.71028251,
                411.66838031, 271.65612259, 246.64015526, 138.28842249,
                243.2849937 , 218.46696429, 378.21048668, 289.04864315,
                161.83926452,  68.06962877, 285.2890697 , 298.05328173,
                178.49149643, 229.81611787, 395.29432446, 383.55011643,
                401.28378546, 141.46165188, 403.53048243, 155.22359805,
                194.71797787, 391.75845968, 377.74151766, 285.9081199 ,
                50.32519898, 326.2868063 , 259.86473339, 360.31036159,
                404.21344621, 299.63852993, 179.32381088,  90.55727011,
                80.91828583,  49.98606852, 199.49816491, 239.50182318,
                340.06342879, 162.27011435, 318.38666731,  77.32577414,
                358.03533754, 200.46274264, 159.66667801, 181.77244339,
                387.13949422, 119.1398619 , 194.35347059, 197.26089076,
                411.39113312,  85.35517451, 235.72717444, 236.52094041,
                157.74786941, 368.2695155 , 193.57571873, 114.73353191,
                228.43695593, 278.53328793, 112.90781775, 176.44027242};
        x[2] = new double[x[0].length];
        Arrays.fill(x[2], 1);

        return x;
    }

    public static ImageExt getOldhouse2A2000000() throws IOException {
        //720 X 480
        String path = ResourceFinder.findTestResourcesDirectory() + sep + MASKSDIR
                + sep + "A2000000.bmp";

        File f = new File(path);
        if (!f.exists()) {
            throw new IOException("could not find file at " + path);
        }

        ImageExt image = ImageIOHelper.readImageExt(path);

        return image;
    }

    public static ImageExt getOldhouse2A2000084() throws IOException {
        //720 X 480
        String path = ResourceFinder.findTestResourcesDirectory() + sep + MASKSDIR
                + sep + "A2000084.bmp";

        File f = new File(path);
        if (!f.exists()) {
            throw new IOException("could not find file at " + path);
        }

        ImageExt image = ImageIOHelper.readImageExt(path);

        return image;
    }
}
