package algorithms.imageProcessing;

import java.awt.FlowLayout;
import java.awt.image.BufferedImage;
import java.awt.image.WritableRaster;
import java.io.IOException;
import java.util.logging.Logger;
import javax.swing.ImageIcon;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.SwingUtilities;

/**
 *
 * @author nichole
 */
public class ImageDisplayer {
       
    protected static Logger log = Logger.getLogger(ImageDisplayer.class.getName());
    
    public static DisposableJFrame displayImage(final String title, final Image img) 
        throws IOException {
        
        if (img == null) {
            return null;
        }
        
        int w = img.getWidth();
        int h = img.getHeight();
        
        final BufferedImage image = new BufferedImage(w, h, 
            BufferedImage.TYPE_INT_RGB);
        
        for (int i = 0; i < w; i++) {
            for (int j = 0; j < h; j++) {
                int rgbValue = img.getRGB(i, j);                
                image.setRGB(i, j, rgbValue);
            }
        }
       
        return displayDisposableImage(title, image);
    }
    
    public static DisposableJFrame displayImage(final String title, 
        final GreyscaleImage img) throws IOException {
        
        if (img == null) {
            return null;
        }
        
        int w = img.getWidth();
        int h = img.getHeight();
        
        final BufferedImage image = new BufferedImage(w, h, 
            BufferedImage.TYPE_BYTE_GRAY);
        
        WritableRaster raster = image.getRaster();
        
        for (int i = 0; i < w; i++) {
            for (int j = 0; j < h; j++) {
                int value = img.getValue(i, j);         
                raster.setSample(i, j, 0, value);
            }
        }
        
        return displayDisposableImage(title, image);
    }

    public static DisposableJFrame displayDisposableImage(final String title, 
        final BufferedImage image) throws IOException {
        
        if (image == null) {
            return null;
        }
       
        try {
            
            boolean isDispatchThread = SwingUtilities.isEventDispatchThread();
              
            log.fine("isDispatchThread=" + isDispatchThread);

            DisposableJFrame rjf = new DisposableJFrame(image, title);
            
            javax.swing.SwingUtilities.invokeAndWait(rjf);

            return rjf;
            
        } catch(Throwable t) {
            t.printStackTrace();
            String msg = t.getMessage();
            log.severe(msg);
        }
        
        return null;
    }
   
    public static class DisposableJFrame implements Runnable {

        protected final JFrame frame;
        protected final BufferedImage image;
        
        public DisposableJFrame(BufferedImage theImage, String title) {
            image = theImage;
            frame = new JFrame(title);
        }
        
        public void dispose() {
            frame.dispose();
        }
        
        @Override
        public void run() {

            ImageIcon icon = new ImageIcon(image);

            JLabel lbl = new JLabel();
            lbl.setIcon(icon);

            frame.setLayout(new FlowLayout());
            frame.setSize(image.getWidth(), image.getHeight());
            frame.add(lbl);
            frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
            frame.setVisible(true);
            frame.toFront();
        }
    }
}
