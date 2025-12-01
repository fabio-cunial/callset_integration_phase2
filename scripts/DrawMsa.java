import java.io.*;
import java.util.*;
import java.awt.*;
import java.awt.image.*;
import java.awt.geom.*;
import javax.imageio.*;
import java.text.*;


/**
 * 
 */
public class DrawMsa {
    
    /**
     * @param 
     */
    public static void main(String[] args) throws IOException {
        final String MSA_TXT = args[0];
        final int N_ROWS = Integer.parseInt(args[1]);
        final int N_COLUMNS = Integer.parseInt(args[2]);
        final String OUTPUT_PNG = args[3];
        
        int x, y;
        int strLength, color;
        String str;
        BufferedReader br;
        BufferedImage image;
        
        image = new BufferedImage(N_COLUMNS,N_ROWS,BufferedImage.TYPE_INT_RGB);
        br = new BufferedReader(new FileReader(MSA_TXT));
        str=br.readLine(); y=-1;
        while (str!=null) {
            if (str.charAt(0)=='>') { y++; str=br.readLine(); continue; }
            strLength=str.length();
            for (x=0; x<str.length(); x++) {
                color=0x00FFFFFF;
                switch (str.charAt(x)) {
                    case 'A': color=0x0058A644; break;
                    case 'C': color=0x000021F5; break;
                    case 'G': color=0x00C77B30; break;
                    case 'T': color=0x00EA3223; break;
                }
                image.setRGB(x,y,color);
            }
            str=br.readLine();
        }
        br.close();
        ImageIO.write(image,"png",new File(OUTPUT_PNG));
    }
    
}