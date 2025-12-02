import java.util.zip.GZIPInputStream;
import java.io.*;


/**
 * Remark: the program handles multiallelic sites, which are common in TRGT.
 */
public class DeNovoRateTrgt {
    
    /**
     * @param args all VCFs are assumed to contain the same set of loci in the
     * same order.
     */
    public static void main(String[] args) throws IOException {
        final String CHILD_VCF_GZ = args[0];
        final String PARENT1_VCF_GZ = args[1];
        final String PARENT2_VCF_GZ = args[2];
        
        final int MIN_SVLEN_1 = 20;
        final int MIN_SVLEN_2 = 50;
        final int QUANTUM = 10000;  // Arbitrary
        
        boolean found;
        int i, j, p;
        int altChildLast, altParent1Last, altParent2Last, numerator, denominator, tmpRow, nRecords;
        String strChild, str1, str2, chr, pos, ref, alt;
        BufferedReader brChild, br1, br2;
        int[][] tmpMatrix, output0, output1, output2;
        String[] tokens, altChild, altParent1, altParent2;
        
        altChild = new String[2]; altParent1 = new String[2]; altParent2 = new String[2];
        altChildLast=-1; altParent1Last=1; altParent2Last=-1;
        brChild = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(CHILD_VCF_GZ))));
        br1 = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(PARENT1_VCF_GZ))));
        br2 = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(PARENT2_VCF_GZ))));
        strChild=brChild.readLine();
        while (strChild.charAt(0)=='#') strChild=brChild.readLine();
        str1=br1.readLine();
        while (str1.charAt(0)=='#') str1=br1.readLine();
        str2=br2.readLine();
        while (str2.charAt(0)=='#') str2=br2.readLine();
        output0 = new int[1][2]; output1 = new int[2][2]; output2 = new int[2][2]; nRecords=0;
        while (strChild!=null) {
            // Collecting the ALTs of the child
            tokens=strChild.split("\t"); chr=tokens[0]; pos=tokens[1]; ref=tokens[3]; alt=tokens[4];
            if (alt.equals(".")) { 
                // Not ALT in the child: not counted.
                strChild=brChild.readLine(); str1=br1.readLine(); str2=br2.readLine(); 
                continue;
            }
            p=alt.indexOf(",");
            if (p!=-1) {
                altChild[0]=alt.substring(0,p);
                altChild[1]=alt.substring(p+1);
                altChildLast=1;
            }
            else {
                altChild[0]=alt;
                altChildLast=0;
            }
            
            // Collecting the ALTs of parent1
            tokens=str1.split("\t"); 
            if (!tokens[0].equals(chr) || !tokens[1].equals(pos) || !tokens[3].equals(ref)) {
                System.err.println("ERROR: the VCFs do not contain the same set of records in the same order.");
                System.exit(1);
            }
            alt=tokens[4];
            if (alt.equals(".")) altParent1Last=-1;
            else {
                p=alt.indexOf(",");
                if (p!=-1) {
                    altParent1[0]=alt.substring(0,p);
                    altParent1[1]=alt.substring(p+1);
                    altParent1Last=1;
                }
                else {
                    altParent1[0]=alt;
                    altParent1Last=0;
                }
            }
            
            // Collecting the ALTs of parent2
            tokens=str2.split("\t");
            if (!tokens[0].equals(chr) || !tokens[1].equals(pos) || !tokens[3].equals(ref)) {
                System.err.println("ERROR: the VCFs do not contain the same set of records in the same order.");
                System.exit(1);
            }
            alt=tokens[4];
            if (alt.equals(".")) altParent2Last=-1;
            else {
                p=alt.indexOf(",");
                if (p!=-1) {
                    altParent2[0]=alt.substring(0,p);
                    altParent2[1]=alt.substring(p+1);
                    altParent2Last=1;
                }
                else {
                    altParent2[0]=alt;
                    altParent2Last=0;
                }
            }
            
            // Counting new ALTs in the child
            for (i=0; i<=altChildLast; i++) {
                // All records
                output0[0][1]++;  // Denominator
                found=false;
                for (j=0; j<=altParent1Last; j++) {
                    if (altParent1[j].equalsIgnoreCase(altChild[i])) found=true;
                }
                for (j=0; j<=altParent2Last; j++) {
                    if (altParent2[j].equalsIgnoreCase(altChild[i])) found=true;
                }
                if (!found) output0[0][0]++;  // Numerator
                tmpMatrix=null; tmpRow=-1;
                // SVs defined by max(|ref|,|alt|)
                if (ref.length()>=MIN_SVLEN_1 || altChild[i].length()>=MIN_SVLEN_1) { tmpMatrix=output1; tmpRow=0; }
                if (ref.length()>=MIN_SVLEN_2 || altChild[i].length()>=MIN_SVLEN_2) { tmpMatrix=output1; tmpRow=1; }
                // SVs defined by abs(|alt|-|ref|)
                if (altChild[i].length()-ref.length()>=MIN_SVLEN_1 || ref.length()-altChild[i].length()>=MIN_SVLEN_1) { tmpMatrix=output2; tmpRow=0; }
                if (altChild[i].length()-ref.length()>=MIN_SVLEN_2 || ref.length()-altChild[i].length()>=MIN_SVLEN_2) { tmpMatrix=output2; tmpRow=1; }
                if (tmpMatrix!=null && tmpRow!=-1) {
                    tmpMatrix[tmpRow][1]++;  // Denominator
                    found=false;
                    for (j=0; j<=altParent1Last; j++) {
                        if (altParent1[j].equalsIgnoreCase(altChild[i])) found=true;
                    }
                    for (j=0; j<=altParent2Last; j++) {
                        if (altParent2[j].equalsIgnoreCase(altChild[i])) found=true;
                    }
                    if (!found) tmpMatrix[tmpRow][0]++;  // Numerator
                }
            }
            
            // Next iteration
            nRecords++;
            if (nRecords%QUANTUM==0) System.err.println("Processed "+nRecords+" records...");
            strChild=brChild.readLine(); str1=br1.readLine(); str2=br2.readLine();
        }
        brChild.close(); br1.close(); br2.close();
        
        // Outputting
        System.out.println("SV definition 1 (max):");
        System.out.println(">="+MIN_SVLEN_1+": "+output1[0][0]+","+output1[0][1]+","+(100.0*output1[0][0])/output1[0][1]+"%");
        System.out.println(">="+MIN_SVLEN_2+": "+output1[1][0]+","+output1[1][1]+","+(100.0*output1[1][0])/output1[1][1]+"%");
        System.out.println("SV definition 2 (diff):");
        System.out.println(">="+MIN_SVLEN_1+": "+output2[0][0]+","+output2[0][1]+","+(100.0*output2[0][0])/output2[0][1]+"%");
        System.out.println(">="+MIN_SVLEN_2+": "+output2[1][0]+","+output2[1][1]+","+(100.0*output2[1][0])/output2[1][1]+"%");
        System.out.println("All records:");
        System.out.println(output0[0][0]+","+output0[0][1]+","+(100.0*output0[0][0])/output0[0][1]+"%");
    }
    
}