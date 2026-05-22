import java.util.*;
import java.util.zip.GZIPInputStream;
import java.io.*;


/**
 * Given a VCF that contains only INS records, the program separates records
 * that are likely DUP based on the output of `UltralongDepthGetBreakpoints`.
 * 
 * Remark: every record in the output DUP VCF has tag INSDUP that indicates that
 * it is the result of INS conversion, and a field INS_ALT that contains the ALT
 * allele of the original INS record.
 * 
 * Remark: the output INS VCF is sorted, but the output DUP VCF is not 
 * necessarily sorted.
 */
public class UltralongInsExtractDups {
    
    /**
     * @param args 4: fixed QUAL value to assign to all INS->DUP records.
     */
    public static void main(String[] args) throws IOException {
        final String INPUT_VCF_GZ = args[0];
        final String BREAKPOINTS_DIR = args[1];
        final String OUTPUT_VCF_INS = args[2];
        final String OUTPUT_VCF_DUP = args[3];
        final int INSDUP_QUAL = Integer.parseInt(args[4]);
        
        int i, p;
        int pos, newPos, newEnd, newLength, nRecords, nDups;
        String str, strPrime, id, info, alt;
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
                if (str.startsWith("#CHROM")) {
                    bwDup.write("##INFO=<ID=INSDUP,Number=0,Type=Flag,Description=\"The record is the result of an INS->DUP conversion\">\n");
                    bwDup.write("##INFO=<ID=INS_ALT,Number=1,Type=String,Description=\"The ALT allele of the original INS record\">\n");
                }
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
                alt=tokens[4]; tokens[4]="<DUP>";
                tokens[5]=INSDUP_QUAL+"";
                info=addOrReplaceInfoField(info,"SVLEN",String.valueOf(newLength));
                info=addOrReplaceInfoField(info,"SVTYPE","DUP");
                info=addOrReplaceInfoField(info,"END",newEnd+"");
                info+=";INSDUP;INS_ALT="+alt;
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