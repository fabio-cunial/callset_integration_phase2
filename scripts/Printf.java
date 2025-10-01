
/**
 * 
 */
public class Printf {
    
    /**
     * @param
     */
    public static void main(String[] args) {
        final int MIN_DEPTH = Integer.parseInt(args[0]);
        final int MAX_DEPTH = Integer.parseInt(args[1]);
        final int MIN_ALT_READS = Integer.parseInt(args[2]);
        
        System.out.println("'GT=\"alt\" && (DP < "+MIN_DEPTH+" || DP > "+MAX_DEPTH+" || AD[*:1] < "+MIN_ALT_READS+")'");
    }
    
}
