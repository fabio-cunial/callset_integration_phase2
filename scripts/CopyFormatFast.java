import java.io.*;


/**
 * 
 */
public class CopyFormatFast {
    
    /**
     * 
     */
    public static void main(String[] args) throws IOException {
        final String CHUNK_FROM = args[0];
        final String CHUNK_TO = args[1];
        final String OUT_CHUNK = args[2];
        final String OUT_HWE = args[3];
        
        final char FIELD_SEPARATOR = '\t';
        final char GT_SEPARATOR = ':';
        final int QUANTUM = 100;  // Arbitrary
        
        boolean onLeft, onRight;
        char c;
        int i, p;
        int lengthFrom, lengthTo, nLines, startFrom, startTo, endFrom, endTo;
        int gt00, gt01, gt11;
        long ms;
        String strFrom, strTo;
        BufferedReader brFrom, brTo;
        BufferedWriter bw, bwHwe;
        
        brFrom = new BufferedReader(new FileReader(CHUNK_FROM));
        brTo = new BufferedReader(new FileReader(CHUNK_TO));
        bw = new BufferedWriter(new FileWriter(OUT_CHUNK));
        bwHwe = new BufferedWriter(new FileWriter(OUT_HWE));
        strFrom=brFrom.readLine(); strTo=brTo.readLine(); nLines=0;
        ms=System.currentTimeMillis();
        while (strTo!=null) {
            gt00=0; gt01=0; gt11=0;
            lengthFrom=strFrom.length(); lengthTo=strTo.length();
            startFrom=0; startTo=0;
            for (i=0; i<=8; i++) startTo=strTo.indexOf(FIELD_SEPARATOR,startTo)+1;
            for (i=0; i<startTo-1; i++) bw.write(strTo.charAt(i));
            do {
                bw.write(FIELD_SEPARATOR);
                // To
                endTo=strTo.indexOf(FIELD_SEPARATOR,startTo);
                if (endTo<0) endTo=lengthTo;
                p=startTo;
                for (i=0; i<9; i++) {
                    p=strTo.indexOf(GT_SEPARATOR,p);
                    if (p>=endTo || p<0) { p=endTo; break; }
                    else p++;
                }
                c='_';
                for (i=startTo; i<p; i++) { c=strTo.charAt(i); bw.write(c); }
                if (p==endTo && c!=GT_SEPARATOR) bw.write(GT_SEPARATOR);
                // From
                endFrom=strFrom.indexOf(FIELD_SEPARATOR,startFrom);
                if (endFrom<0) endFrom=lengthFrom;
                p=startFrom;
                for (i=0; i<9; i++) p=strFrom.indexOf(GT_SEPARATOR,p)+1;
                for (i=p; i<endFrom; i++) bw.write(strFrom.charAt(i));
                // Updating genotype counts using $strTo$.
                p=strTo.indexOf(GT_SEPARATOR,startTo);
                if (p==startTo+3) {
                    onLeft=strTo.charAt(startTo)=='1';
                    onRight=strTo.charAt(startTo+2)=='1';
                    if (onLeft==onRight) {
                        if (onLeft) gt11++;
                        else gt00++;
                    }
                    else gt01++;
                }
                // Next sample
                startFrom=endFrom+1; startTo=endTo+1;
            } while (startFrom<lengthFrom && startTo<lengthTo);
            if ((startFrom<lengthFrom)!=(startTo<lengthTo)) {
                System.err.println("ERROR: the two lines do not contain the same number of samples.");
                System.err.println("strFrom: "+strFrom);
                System.err.println("strTo: "+strTo);
                System.exit(1);
            }
            bw.newLine(); 
            if (gt00+gt01+gt11>0) bwHwe.write(gt00+","+gt01+","+gt11+"\n");
            nLines++;
            if (nLines%QUANTUM==0) System.err.println("Processed "+nLines+" lines in "+((double)(System.currentTimeMillis()-ms)/1000)+"s");
            strFrom=brFrom.readLine(); strTo=brTo.readLine();
            if ((strFrom==null)!=(strTo==null)) {
                System.err.println("ERROR: the two files do not contain the same number of lines.");
                System.exit(1);
            }
        }
        brFrom.close(); brTo.close(); bw.close(); bwHwe.close();
    }
    
}