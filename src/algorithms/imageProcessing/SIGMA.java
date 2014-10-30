package algorithms.imageProcessing;

/**
 *
 * @author nichole
 */
public enum SIGMA {
    
    //for result where convolution is mult by 1, sigma=0.42466090014400953
    // which is close to 0.5
    /*ZEROPOINTFOURTWOSEVEN,*/ 
    ZEROPOINTFIVE, ZEROPOINTSEVENONE,
    ONE, ONESQRT2, ONEPOINTFIVE,
    TWO, ONEPOINTFIVESQRT2, TWOSQRT2,
    THREE,
    FOUR, FOURSQRT2,
    EIGHT, EIGHTSQRT2,
    SIXTEEN, SIXTEENSQRT2,
    THIRTYTWO, THIRTYTWOSQRT2,
    SIXTYFOUR, SIXTYFOURSQRT2,
    ONEHUNDREDANDTWENTYEIGHT, ONEHUNDREDANDTWENTYEIGHTSQRT2,
    TWOHUNDREDANDFIFTYSIX
    ;
    
    public static SIGMA resolve(String name) {
        for (SIGMA vs : SIGMA.values()) {
            if (vs.name().equals(name)) {
                return vs;
            }
        }
        return null;
    }
    
    public static float getValue(SIGMA sigma) {
    
        /*if (sigma.ordinal() == SIGMA.ZEROPOINTFOURTWOSEVEN.ordinal()) {
            return 0.42466090014400953f;
        } else*/ if (sigma.ordinal() == SIGMA.ZEROPOINTFIVE.ordinal()) {
            return 0.5f;
        } else if (sigma.ordinal() == SIGMA.ZEROPOINTSEVENONE.ordinal()) {
            return (float) (1./Math.sqrt(2));
        } else if (sigma.ordinal() == SIGMA.ONE.ordinal()) {
            return 1;
        } else if (sigma.ordinal() == SIGMA.ONESQRT2.ordinal()) {
            return (float)Math.sqrt(2);
        } else if (sigma.ordinal() == SIGMA.ONEPOINTFIVE.ordinal()) {
            return 1.5f;
        } else if (sigma.ordinal() == SIGMA.TWO.ordinal()) {
            return 2;
        } else if (sigma.ordinal() == SIGMA.ONEPOINTFIVESQRT2.ordinal()) {
            return 1.5f*(float)Math.sqrt(2);
        } else if (sigma.ordinal() == SIGMA.TWOSQRT2.ordinal()) {
            return 2*(float)Math.sqrt(2);
        } else if (sigma.ordinal() == SIGMA.THREE.ordinal()) {
            return 3;
        } else if (sigma.ordinal() == SIGMA.FOUR.ordinal()) {
            return 4;
        } else if (sigma.ordinal() == SIGMA.FOURSQRT2.ordinal()) {
            return 4*(float)Math.sqrt(2);
        } else if (sigma.ordinal() == SIGMA.EIGHT.ordinal()) {
            return 8;
        } else if (sigma.ordinal() == SIGMA.EIGHTSQRT2.ordinal()) {
            return 8*(float)Math.sqrt(2);
        } else if (sigma.ordinal() == SIGMA.SIXTEEN.ordinal()) {
            return 16;
        } else if (sigma.ordinal() == SIGMA.SIXTEENSQRT2.ordinal()) {
            return 16*(float)Math.sqrt(2);
        } else if (sigma.ordinal() == SIGMA.THIRTYTWO.ordinal()) {
            return 32;
        } else if (sigma.ordinal() == SIGMA.THIRTYTWOSQRT2.ordinal()) {
            return 32*(float)Math.sqrt(2);
        } else if (sigma.ordinal() == SIGMA.SIXTYFOUR.ordinal()) {
            return 64;
        } else if (sigma.ordinal() == SIGMA.SIXTYFOURSQRT2.ordinal()) {
            return 64*(float)Math.sqrt(2);
        } else if (sigma.ordinal() == 
            SIGMA.ONEHUNDREDANDTWENTYEIGHT.ordinal()) {
            return 128;
        } else if (sigma.ordinal() == 
            SIGMA.ONEHUNDREDANDTWENTYEIGHTSQRT2.ordinal()) {
            return 128*(float)Math.sqrt(2);
        } else if (sigma.ordinal() == 
            SIGMA.TWOHUNDREDANDFIFTYSIX.ordinal()) {
            return 256;
        }
        
        throw new IllegalArgumentException("Programmer Error!  Didn't include " 
            + sigma.toString());
    }
    
    public static SIGMA increment(SIGMA sigma) {
        
        if (sigma.ordinal() == (SIGMA.values().length - 1)) {
            return null;
        }
        
        return SIGMA.values()[sigma.ordinal() + 1];
    }
  
    public static SIGMA divideBySQRT2(SIGMA sigma) {
        
        if (sigma.ordinal() == SIGMA.FOUR.ordinal()) {
            return SIGMA.TWOSQRT2;
        } else if (sigma.ordinal() == SIGMA.TWOSQRT2.ordinal()) {
            return SIGMA.TWO;
        } else if (sigma.ordinal() == SIGMA.TWO.ordinal()) {
            return SIGMA.ONESQRT2;
        } else if (sigma.ordinal() == SIGMA.ONESQRT2.ordinal()) {
            return SIGMA.ONE;
        } else if (sigma.ordinal() == SIGMA.ONE.ordinal()) {
            return SIGMA.ZEROPOINTSEVENONE;
        } else if (sigma.ordinal() == SIGMA.ONEPOINTFIVE.ordinal()) {
            return null;
            
        } else if (sigma.ordinal() == SIGMA.THREE.ordinal()) {
            return SIGMA.ONEPOINTFIVESQRT2;
        } else if (sigma.ordinal() == SIGMA.ONEPOINTFIVESQRT2.ordinal()) {
            return SIGMA.ONEPOINTFIVE;
        } else if (sigma.ordinal() == SIGMA.ZEROPOINTSEVENONE.ordinal()) {
            return SIGMA.ZEROPOINTFIVE;
        } else if (sigma.ordinal() == SIGMA.ZEROPOINTFIVE.ordinal()) {
            return null;
        }
   
        throw new IllegalArgumentException("have not written the division for "
            + sigma.toString());
    }
    
    public static SIGMA increaseToFactorBySQRT2(float sigma) {
        
        /*if (Math.abs(sigma - 0.4247f) < 0.01) {
            //0.42466090014400953f
            return SIGMA.ZEROPOINTFOURTWOSEVEN;
        } else*/ if (Math.abs(sigma - 0.5f) < 0.01f) {
            return SIGMA.ZEROPOINTFIVE;
        } else if (Math.abs(sigma - 0.707106f) < 0.01f) {
            return SIGMA.ZEROPOINTSEVENONE;
        } else if (sigma == 1) {
            return SIGMA.ONE;
        } else if (sigma == ((float)Math.sqrt(2))) {
            return SIGMA.ONESQRT2;
        } else if (sigma == 2) {
            return SIGMA.TWO;
        } else if (sigma == (2*(float)Math.sqrt(2))) {
            return SIGMA.TWOSQRT2;
        } else if (sigma == 4) {
            return SIGMA.FOUR;
        } else if (sigma == (4*(float)Math.sqrt(2))) {
            return SIGMA.FOURSQRT2;
        } else if (sigma == 8) {
            return SIGMA.EIGHT;
        } else if (sigma == (8*(float)Math.sqrt(2))) {
            return SIGMA.EIGHTSQRT2;
        } else if (sigma == 16) {
            return SIGMA.SIXTEEN;
        } else if (sigma == (16*(float)Math.sqrt(2))) {
            return SIGMA.SIXTEENSQRT2;
        } else if (sigma == 32) {
            return SIGMA.THIRTYTWO;
        } else if (sigma == (32*(float)Math.sqrt(2))) {
            return SIGMA.THIRTYTWOSQRT2;
        } else if (sigma == 64) {
            return SIGMA.SIXTYFOUR;
        } else if (sigma == (64*(float)Math.sqrt(2))) {
            return SIGMA.SIXTYFOURSQRT2;
        } else if (sigma == 128) {
            return SIGMA.ONEHUNDREDANDTWENTYEIGHT;
        } else if (sigma == (128*(float)Math.sqrt(2))) {
            return SIGMA.ONEHUNDREDANDTWENTYEIGHTSQRT2;
        } else if (sigma == 256) {
            return SIGMA.TWOHUNDREDANDFIFTYSIX;
        } else if (sigma > 256) {
            return null;
        
        } else if (Math.abs(sigma - 0.5) < 0.01) {
            return SIGMA.ZEROPOINTFIVE;
        } else if (Math.abs(sigma - 1.5f*Math.sqrt(2)) < 0.01) {
            return SIGMA.ONEPOINTFIVESQRT2;
        } else if (sigma == 3) {
            return SIGMA.THREE;
        }
        
        throw new IllegalArgumentException("have not written the increment for "
            + sigma);
    }
    
}
