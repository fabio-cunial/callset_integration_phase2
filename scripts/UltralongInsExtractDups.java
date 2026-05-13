import java.util.*;
import java.util.zip.GZIPInputStream;
import java.io.*;


/**
 * Given a VCF that contains only INS records, the program separates records
 * that are likely DUP based on the output of `UltralongDepthGetBreakpoints`.
 * 
 * Remark: the output INS VCF is sorted, but the output DUP VCF is not 
 * necessarily sorted.
 */
public class UltralongInsExtractDups {
    
    /**
     * @param args
     */
    public static void main(String[] args) throws IOException {
        final String INPUT_VCF_GZ = args[0];
        final String BREAKPOINTS_DIR = args[1];
        final String OUTPUT_VCF_INS = args[2];
        final String OUTPUT_VCF_DUP = args[3];
        
        int i, p;
        int pos, newPos, newEnd, newLength, nRecords, nDups;
        String str, strPrime, id, info;
        BufferedReader br, brPrime;
        BufferedWriter bwIns, bwDup;
        String[] tokens;

        bwIns = new BufferedWriter(new FileWriter(OUTPUT_VCF_INS));
        bwDup = new BufferedWriter(new FileWriter(OUTPUT_VCF_DUP));
        br = new BufferedReader( new InputStreamReader( (INPUT_VCF_GZ.length()>=7&&INPUT_VCF_GZ.substring(INPUT_VCF_GZ.length()-7).equalsIgnoreCase(".vcf.gz")) ? new GZIPInputStream(new FileInputStream(INPUT_VCF_GZ)) : new FileInputStream(INPUT_VCF_GZ) ) );
        str=br.readLine(); nRecords=0; nDups=0;
        while (str!=null) {
            if (str.charAt(0)=='#') {
                bwIns.write(str+"\n");
                bwDup.write(str+"\n");
                str=br.readLine();
                continue;
            }
            nRecords++;
            tokens=str.split("\t");
            id=tokens[2];
            info=tokens[7];
            brPrime = new BufferedReader(new FileReader(BREAKPOINTS_DIR+"/"+id+"_breakpoints.tsv"));
            strPrime=brPrime.readLine();
            brPrime.close();
            if (strPrime==null) bwIns.write(str+"\n");
            else {
                nDups++;
                p=strPrime.indexOf("\t");
                newPos=Integer.parseInt(strPrime.substring(0,p));
                newEnd=Integer.parseInt(strPrime.substring(p+1));
                newLength=newEnd-newPos;
                tokens[1]=newPos+"";
                tokens[4]="<DUP>";
                info=addOrReplaceInfoField(info,"SVLEN",String.valueOf(newLength));
                info=addOrReplaceInfoField(info,"SVTYPE","DUP");
                info=addOrReplaceInfoField(info,"END",newEnd+"");
                tokens[7]=info;
                bwDup.write(tokens[0]);
                for (i=1; i<tokens.length; i++) bwDup.write("\t"+tokens[i]);
                bwDup.write("\n");
            }
            str=br.readLine();
        }
        br.close(); bwIns.close(); bwDup.close();
        System.err.println("Created "+nDups+" DUPs out of "+nRecords+" INS records ("+String.format("%.2f",(100.0*nDups)/nRecords)+"%).");
    }


    private static final String addOrReplaceInfoField(String info, String field, String newValue) {
		final int FIELD_LENGTH = field.length()+1;
        int p, q;
        
        if (info.equals(".")) return field+"="+newValue;
        p=-FIELD_LENGTH;
        do { p=info.indexOf(field+"=",p+FIELD_LENGTH); }
        while (p>0 && info.charAt(p-1)!=';');
		if (p<0) return info+";"+field+"="+newValue;
		q=info.indexOf(";",p+FIELD_LENGTH);
        return info.substring(0,p+FIELD_LENGTH)+newValue+(q>=0?info.substring(q):"");
	}

}