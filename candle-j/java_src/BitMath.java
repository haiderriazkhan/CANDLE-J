package bit.math;

public class BitMath {
    public static int lp2(int num) {
        if (num == 0) return 0;
        int result = 1;
        while ((result < num) && (result != 0)) {
            result <<= 1;
        }
        return result;
    }
}
