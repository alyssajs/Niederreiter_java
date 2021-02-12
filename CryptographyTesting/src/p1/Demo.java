package p1;

import java.util.Scanner;

public class Demo 
{
   public static void main(String[] args)
   {
      Scanner scan = new Scanner(System.in);
      
      int fieldExponent =  5;
      int numErrors = 3;
      int supportSize = (int)(Math.pow(2, fieldExponent));
      //int codeLength = 0;
      
      System.out.println("Using parameters:"
            + "\n Field exponent: " + fieldExponent
            + "\n Number of errors to correct: " + numErrors);  
      
     int[] toEncrypt = McEliece.getRandomErrorVector(2,supportSize -  fieldExponent*numErrors);
      
      
 
      System.out.println("Encrypting: ");
      McEliece.print(toEncrypt);
      

   
      int[] encrypted = McEliece.encrypt(toEncrypt, numErrors, fieldExponent, supportSize);
      
      System.out.println("Goppa Polynomial: " + McEliece.irredPoly.toString());
      System.out.println("Field Polynomial: " + McEliece.gF.gf_irredPoly.toString());
      System.out.println("Support size: " + McEliece.support.length);
      System.out.println("Minimum distance >= " + (2 * McEliece.irredPoly.degree + 1));
      System.out.println("Error correction capacity: " + McEliece.irredPoly.degree);
      scan.nextLine();
     
      
      System.out.println("Parity check:");
      McEliece.print(McEliece.parityCheck);
      if(!McEliece.isSystematic(McEliece.parityCheck, McEliece.parityCheck[0].length - McEliece.parityCheck.length))
      {
         scan.nextLine();
      }

      
      System.out.println("Public key: ");
      McEliece.print(McEliece.pubKey);

      System.out.println("Encrypted:");
      McEliece.print(encrypted);
      //DECRYPTION
 

      int[] decrypted = McEliece.decrypt(encrypted, McEliece.parityCheck);
      System.out.println("Message = Decrypted?");
      System.out.println(McEliece.vecEquals(toEncrypt, decrypted));
      
      
      }
      
   
   
}
