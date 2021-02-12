package p1;
import java.math.BigInteger;
import java.util.*; 


public class Tester {

   public static void main(String[] args) 
   {
      Scanner scan = new Scanner(System.in);
      int fieldExponent =  4;
      int numErrors = 2;
      int supportSize = (int)(Math.pow(2, fieldExponent));
      //int codeLength = 0;
      
      System.out.println("Using parameters:"
            + "\n Field exponent: " + fieldExponent
            + "\n Number of errors to correct: " + numErrors);  

     //int[] toEncrypt = McEliece.getRandomErrorVector(1,supportSize -  fieldExponent*numErrors);
     int[] toEncrypt = new int[] { 0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1};
      System.out.println("Original");
     McEliece.print(toEncrypt);
     //int[] enc = McEliece.encrypt(toEncrypt, numErrors, fieldExponent, supportSize);
     
     int[] encoded = Niederreiter.encode(toEncrypt, numErrors, supportSize);
     System.out.println("CW encoded: ");
     McEliece.print(encoded);
     int[] niedEnc = Niederreiter.encrypt(encoded,Niederreiter.createParityCheck(fieldExponent,numErrors, supportSize));
System.out.println("PC:");
McEliece.print(Niederreiter.parityCheck);
     
     System.out.println("Nied encrypt:");
     McEliece.print(niedEnc);

 
     
     
     System.out.println("Decrypting:");
     int[] dec = Niederreiter.decrypt(niedEnc);
     System.out.println("Support:");
     McEliece.print(Niederreiter.support);

   System.out.println("permInv");
   McEliece.print(McEliece.invert(Niederreiter.permutation));

    
     
     
      
    //int[] toEncrypt = McEliece.getRandomErrorVector(numErrors, supportSize);  
 /*
      System.out.println("Encrypting: ");
      McEliece.print(toEncrypt);
      
      int[] toEncode = new int[] {1,0,1,0,1};
      int[] result = Niederreiter.encode(toEncode, numErrors, supportSize);
      System.out.println("encoded: ");
      McEliece.print(result);
      

      scan.nextLine(); 
      //int[] encrypted = McEliece.encrypt(toEncrypt, numErrors, fieldExponent, supportSize);
      int[][] pC = Niederreiter.createParityCheck(fieldExponent, numErrors, supportSize);
      int[] encrypted = Niederreiter.encrypt(result, pC);
      System.out.println("encrypted: ");
      McEliece.print(encrypted);
    
      System.out.println("decoded: ");
      int[] decodedTest = Niederreiter.decodeCW(numErrors, supportSize, Niederreiter.zeroes);
      McEliece.print(decodedTest);
       
      
      int[] dec = Niederreiter.decrypt(encrypted);
      System.out.println("decrypted: ");
      McEliece.print(dec);
      
      /*
      
      GaloisField gf = new GaloisField(3);
      gfPoly[] coeffs = new gfPoly[] { gfPoly.gfOne, new gfPoly("10", 1), 
            new gfPoly("110", 2)};
      
      Poly p1 = new Poly(coeffs);
      
      Poly trace = McEliece.getTracePolynomial();
      
      System.out.println("Factoring: " + p1.toString());
      
      System.out.println("Trace: " + trace.toString());
      List<gfPoly> roots = McEliece.runBerl(p1, 0, trace);
      System.out.print("Roots: ");
      for(int i = 0; i < roots.size(); i++)
      {
         System.out.print(roots.get(i).toString());  
         if(i != roots.size() - 1)
         {
            System.out.print(", ");
         }
      }
      /*
      int[][] toInvert = new int[][] {{1,0,1,0},{0,1,0,1},{1,0,1,1},{0,1,1,0}};
      McEliece.invert(toInvert);
      */
   }
   
   
}
      
      
