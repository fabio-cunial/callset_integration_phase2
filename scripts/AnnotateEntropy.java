import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.zip.GZIPInputStream;
import java.io.*;


/**
 * 
 */
public class AnnotateEntropy {
    
    /**
     * Remark: the output VCF is printed to STDOUT.
     *
     * @param args 1: comma-separated list of entropy orders.
     */
    public static void main(String[] args) throws IOException {
        final String INPUT_VCF_GZ = args[0];
        final String ORDERS = args[1];
        
        final int sigma = 4;
        
        int i;
        int refLength, altLength;
        double entropy;
        String str;
        StringBuilder sequence;
        BufferedReader br;
        int[] orders;
        double[] tmpCounts;
        StringBuilder[] tmpSbs;
        String[] tokens;
        HashMap<String,StringBuilder> tmpMap;
        
        tokens=ORDERS.split(",");
        orders = new int[tokens.length];
        for (i=0; i<tokens.length; i++) orders[i]=Integer.parseInt(tokens[i]);
        sequence = new StringBuilder();
        tmpMap = new HashMap();
        tmpCounts = new double[sigma];
        tmpSbs = new StringBuilder[sigma];
        for (i=0; i<sigma; i++) tmpSbs[i] = new StringBuilder();
        br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(INPUT_VCF_GZ))));
        str=br.readLine();
        while (str!=null) {
            if (str.charAt(0)=='#') {
                if (str.substring(0,6).equalsIgnoreCase("#CHROM")) {
                    for (i=0; i<orders.length; i++) System.out.println("##INFO=<ID=H"+orders[i]+",Number=1,Type=Float,Description=\""+orders[i]+"-order empirical entropy\">");
                }
                System.out.println(str);
                str=br.readLine();
                continue;
            }
            tokens=str.split("\t");
            refLength=tokens[3].length(); altLength=tokens[4].length();
            sequence.delete(0,sequence.length());
            if (refLength>1 && altLength==1) sequence.append(tokens[3]);
            else if (refLength==1 && altLength>1) sequence.append(tokens[4]);
            else if (refLength>1 && altLength>1) {
                // Arbitrary decision
                sequence.append(tokens[4]);
            }
            for (i=0; i<orders.length; i++) {
                entropy=entropy(orders[i],sequence,tmpMap,tmpCounts,tmpSbs);
                tokens[7]=tokens[7]+";H"+orders[i]+"="+entropy;
            }
            System.out.print(tokens[0]);
            for (i=1; i<tokens.length; i++) System.out.print("\t"+tokens[i]);
            System.out.println();
            str=br.readLine();
        }
        br.close();
    }
    
    
    /**
     * Naive, slow implementation.
     *
     * @return -1 if `sequence` is empty or has length `k` or less.
     */
    private static final double entropy(int k, StringBuilder sequence, HashMap<String,StringBuilder> tmpMap, double[] tmpCounts, StringBuilder[] tmpSbs) {
        final int sigma = 4;
        final double n = sequence.length();
        
        char b, c;
        int i;
        double out;
        String kmer;
        StringBuilder sb;
        
        if (n==0) return -1;
        else if (n<=k) return -1;
        else if (k==0) return entropy0(sequence,tmpCounts);
        else if (k==1) {
            for (i=0; i<sigma; i++) tmpSbs[i].delete(0,tmpSbs[i].length());
            for (i=0; i<n-1; i++) {
                b=sequence.charAt(i);
                c=sequence.charAt(i+1);
                if (b=='a' || b=='A') tmpSbs[0].append(c);
                else if (b=='c' || b=='C') tmpSbs[1].append(c);
                else if (b=='g' || b=='G') tmpSbs[2].append(c);
                else if (b=='t' || b=='T') tmpSbs[3].append(c);
            }
            out=0;
            for (i=0; i<sigma; i++) {
                if (tmpSbs[i].length()!=0) out+=tmpSbs[i].length()*entropy0(tmpSbs[i],tmpCounts);
            }
            return out/n;
        }
        else {
            tmpMap.clear();
            for (i=0; i<n-k; i++) {
                kmer=sequence.substring(i,i+k);
                c=sequence.charAt(i+k);
                if (!tmpMap.containsKey(kmer)) {
                    sb = new StringBuilder();
                    sb.append(c);
                    tmpMap.put(kmer,sb);
                }
                else tmpMap.get(kmer).append(c);
            }
            Iterator<Map.Entry<String,StringBuilder>> iterator = tmpMap.entrySet().iterator();
            out=0;
            while (iterator.hasNext()) {
                sb=iterator.next().getValue();
                out+=sb.length()*entropy0(sb,tmpCounts);
            }
            return out/n;
        }
    }
    
    
    
    // -----------> Implement also a faster order-1 and order-2
    
    
    
    /**
     * Naive, slow implementation.
     */
    private static final double entropy0(StringBuilder sequence, double[] tmp) {
        final int sigma = 4;
        final double n = sequence.length();
        
        char c;
        int i;
        double out;

        for (i=0; i<sigma; i++) tmp[i]=0;
        for (i=0; i<n; i++) {
            c=sequence.charAt(i);
            if (c=='a' || c=='A') tmp[0]++;
            else if (c=='c' || c=='C') tmp[1]++;
            else if (c=='g' || c=='G') tmp[2]++;
            else if (c=='t' || c=='T') tmp[3]++;
        }
        out=0;
        for (i=0; i<sigma; i++) {
            if (tmp[i]!=0) out+=tmp[i]*Math.log(n/tmp[i]);
        }
        return out/n;
    }
    
}