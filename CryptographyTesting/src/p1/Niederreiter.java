package p1;

import java.util.List;
import java.util.Scanner;

public class Niederreiter 
{
   public static int extDegree;
   
   public static int numErrors;
   
   public static gfPoly[] support;
   
   public static Poly irredPoly;
   
   public static Scanner scan = new Scanner(System.in);
   
   public static GaloisField gF;
   
   public static int[][] parityCheck;
   
   public static int[][] permutation;
   
   public static int[][] scramble;
      
   public static Poly[][] polyParity;
   
   public static int[][] unreducedParity;
   
   public static int[] zeroes;
   
   public static int currMsgLoc;
   
   //source: compact constant weight encoding engines
   static int[] encodeCW(int[] msg, int weight, int length)
   {
      int delta = 0;
      int index = 1;
      int[] zeroes = new int[weight];
      int i;
      currMsgLoc = 0;
      int d;
      
      while(weight > 0)
      {
         if(length <= weight)
         {
            //System.out.println("length <= weight, so in here with " + length + " and weight " + weight);
            weight--;
            length--;
           // System.out.println("and new delta is " + delta);
            zeroes[index - 1] = delta;
            delta = 0;
            index++;
         }
         else
         {
            d = calcBestD(weight, length);
            //System.out.println("and d is now " + d + " with weight " + weight + " and length " + length);
           // System.out.println("\n getting entry at " + currMsgLoc + "\n");

            if(msg[currMsgLoc] == 1) 
            {
               currMsgLoc++;
               length -= d;
               delta += d;
            }
            else
            {
               currMsgLoc++;
               i = decodefd(msg, d);
               //System.out.println("Got " + i + " from decodefd");
               zeroes[index - 1] = delta + i;
               length -= (i+1);
               weight--;
               delta = 0;
               index++;
            }
         }
      }
      
      for(int test = 0; test < zeroes.length; test++)
      {
         System.out.print(zeroes[test] + " ");
      }
      System.out.println();
      return zeroes;
   }
   
   static int calcBestD(int weight, int length)
   {
      double term1 = length - (double)(weight-1)/2;
      double pow = (double)1/weight;

      double term2 = 1 - (double)1/Math.pow(2,pow);

      double logVal = Math.log(term1*term2)/Math.log(2);
      int d = (int)Math.pow(2, Math.ceil(logVal));
      return d;
   }
   
   static int[] decodeCW(int weight, int length, int[] zeroes)
   {
      System.out.println("came in with weight: " + weight + " and length " + length);
      int index = 0;
      int[] message = new int[1];
      int d;
      while(weight != 0 && length > weight)
      {
         d = calcBestD(weight, length);
         if(zeroes[index] >= d)
         {
            length -= d;
            zeroes[index] -= d;
            message = addAndResize(message, 1);
            System.out.println("added one thing to message");

         }
         else
         {
            message = addAndResize(message, 0);
            System.out.println("added 0 to message");
            int[] newEntries = encodefd(zeroes[index], d);
            message = addAndResize(message, newEntries);
            System.out.println("added " + newEntries.length + " things to message");
            length -= (zeroes[index] + 1);
            weight--;
            index++;
            
         }
      }
      
      return message;
   }
   
   
   static int[] addAndResize(int[] orig, int newEntry)
   {
      int[] newArray = new int[orig.length + 1];
      for(int loc = 0; loc < orig.length; loc++)
      {
         newArray[loc] = orig[loc];
      }
      newArray[newArray.length - 1] = newEntry;
      return newArray;
      
   }
   
   static int[] addAndResize(int[] orig, int[] newEntries)
   {
      int[] newArray = new int[orig.length + newEntries.length];
      for(int loc = 0; loc < newArray.length; loc++)
      {
         if(loc < orig.length)
         {
            newArray[loc] = orig[loc];
         } 
         else
         {
            newArray[loc] = newEntries[loc - orig.length];
         }
      }
     
      return newArray;
   }
   
   //source: sendrier cw encoding
   static int[] encodefd(int delta, int d)
   {
      int u =(int)Math.ceil(Math.log(d)/Math.log(2));
      int check = (int)(Math.pow(2, u) - d);
      if(delta < check)
      {
         u--;
      }
      else
      {
         delta += check;
      }
      
      //return the u least significant bits of delta written base 2
      int[] result = new int[u];
      for(int digit = 0; digit < u; digit++)
      {
         result[u - digit - 1] = delta % d;
         delta = delta / 2;
      }
      for(int t = 0; t < result.length; t++)
      {
         System.out.print(result[t] + " ");
      }
      System.out.println();
      
      return result;
      
   }
   
   static int decodefd(int[] msg, int d)
   {
      int u = (int)Math.ceil(Math.log(d)/Math.log(2));
      int delta = 0;
      /*
      for(int loc = u-1; loc > 0; loc--)
      {
         delta += Math.pow(2, loc) * msg[u - 1 - loc];
      }
      */
      int tempMsgLoc = currMsgLoc;
      int count = 0;
      for(int loc = tempMsgLoc; loc < tempMsgLoc + u - 1; loc++)
      {
    
         delta += Math.pow(2, u - count - 2) * msg[loc];
         currMsgLoc++;
         count++; 
      }
      //System.out.println("Got " + delta + " from delta");
      if(delta >= (Math.pow(2, u) - d))
      {
         
         delta = 2 * delta + msg[currMsgLoc] - (int)Math.pow(2, u) + d;
         currMsgLoc++;
      }
      return delta;
            
      
   }

   
   static int[] encode(int[] msg, int weight, int length)
   {
      zeroes = encodeCW(msg, weight, length);
      int[] encoded = new int[length];
      int current1 = 0;
      for(int i = 0; i < weight; i++)
      {
         current1 = current1 + zeroes[i];
         if(i != 0)
         {
            current1++;
         }
         encoded[current1] = 1;
     
      }
      
      return encoded; 
   }
   
   public static int[] encrypt(int[] message, int[][] parityCheck)
   {
   
      //placeholder matrices to calculate invertible matrix
      int dim = parityCheck.length;
      int[][] mA = new int[dim][dim];
      int[][] mT = new int[dim][dim];
     System.out.println("Came in with PC dims " + parityCheck.length + " x " + parityCheck[0].length);
      
      //calculate random invertible matrix 
      scramble = McEliece.getRandMatrix(mA, mT, dim);
     
      System.out.println("scramble dims: " + scramble.length + " x " + scramble[0].length);
      
      //calculate permutation matrix 
      permutation = McEliece.createPermMatrix(parityCheck[0].length);
      System.out.println("permutation dims: " + permutation.length + " x " + permutation[0].length);
      
      int[][] messageMatTransp = new int[message.length][1];
      
      
      for(int matRow = 0; matRow < message.length; matRow++)
      {
         messageMatTransp[matRow][0] = message[matRow];
      }
      
      //multiply S * parityCheck * permutation * messageTranspose
      int[][] encryptedMat = McEliece.multiply(scramble, parityCheck);

      encryptedMat = McEliece.multiply(encryptedMat, permutation);
System.out.println("Public key:");
McEliece.print(encryptedMat);
      encryptedMat = McEliece.multiply(encryptedMat, messageMatTransp);
     System.out.println("messageMatTransp dims " + messageMatTransp.length + " x " + messageMatTransp[0].length);
      
      int[] enc = new int[encryptedMat.length];
      for(int loc = 0; loc < enc.length; loc++)
      {
         enc[loc] = encryptedMat[loc][0];
      } 
      
      return enc;
   }
   
   public static int factorial(int num)
   {
      if(num <= 2)
      {
         return num;
      }
      return num * factorial(num - 1);
   }
   
   static int[][] createParityCheck(int extDegreeIn, int numErrorsIn, int suppSize)
   {
      System.out.println("Creating pc with numErrors " + numErrorsIn + " and extDegreeIn " + extDegreeIn);
      parityCheck = McEliece.getParityCheck(extDegreeIn, numErrorsIn, suppSize);
      irredPoly = McEliece.irredPoly;
      extDegree = extDegreeIn;
      numErrors = numErrorsIn;
      support = McEliece.support;
      gF = McEliece.gF;
      polyParity = McEliece.polyParity;  
      unreducedParity = McEliece.unreducedParity;
      return parityCheck;
   }
   
   static int[] decrypt(int[] ciphertext)
   {
      //multiply by random matrix inverse
      int[][] randInverse = McEliece.invert(scramble);
      // int[][] decryptedTemp = new int[1][ciphertext.length];
      
      
      //decryptedTemp[0] = ciphertext;
      
     // decryptedTemp = McEliece.multiply(decryptedTemp, randInverse);
      
      int[][] decryptedTemp = vecTranspose(ciphertext);
      decryptedTemp = McEliece.multiply(randInverse, decryptedTemp);
     McEliece.print(decryptedTemp);
   
      
      int[][] decryptedNext = new int[1][parityCheck[0].length];
   
      //use syndrome decoder (patterson/BTA??)
      
      int[] unscrambled = new int[decryptedTemp.length];
      for(int loc = 0; loc < unscrambled.length; loc++)
      {
         unscrambled[loc] = decryptedTemp[loc][0];
      }
      decryptedNext[0] = pattersonDecode(unscrambled, parityCheck);
      
      //then multiply by permInverse
      int[][] permInverse = McEliece.invert(permutation);
      
      decryptedNext = McEliece.multiply(decryptedNext, permInverse);
      
      return decryptedNext[0];
   }
   

   public static int[] pattersonDecode(int[] codeWord, int[][] parityCheck)
   {
      System.out.println("codeWord.length+ " + codeWord.length);
      System.out.println("support size " + support.length);
      
      System.out.println("in here with: ");
      McEliece.print(codeWord);
      //int[][] scrambleInv = McEliece.invert(scramble);

     
     // int[][] codeTemp = new int[1][codeWord.length];
      //codeTemp[0] = codeWord;
      //System.out.println("codeTemp dims " + codeTemp.length + ", " + codeTemp[0].length);
     // System.out.println("sI dims: " + scrambleInv.length + ", " + scrambleInv[0].length);
      //codeTemp = McEliece.multiply(codeTemp, scrambleInv);

      //codeWord = codeTemp[0];
      
      int[] decoded = new int[parityCheck[0].length];
   
      //Multiply ciphertext by parity check matrix to get syndrome polynomial
      Poly syndrome = getSyndrome(codeWord, numErrors, extDegree);
      System.out.println("Syndrome: " + syndrome.toString());      

      //Invert syndrome polynomial mod goppa polynomial 
      Poly synInverse = Poly.getModularInverse(syndrome, irredPoly, gF.gf_irredPoly).get(1);
System.out.println("Really inverse? " + Poly.multiply(syndrome, synInverse, gF.gf_irredPoly, irredPoly));
      //add f(x) = x to inverse
      //Polynomial add = new Polynomial(synInverse.coefficients + 2);
      Poly add = Poly.add(synInverse, Poly.xPoly);
      //calculate sqrt(x + syndrome inverse)
      Poly sqrt = Poly.calcSqrt(add, irredPoly, gF.gf_irredPoly.degree, gF.gf_irredPoly);
      System.out.println("got sqrt: " + sqrt.toString());
      //System.out.println("sqrt is: " + sqrt.toString());
      //System.out.println("before reducing: " + Poly.toPower(sqrt, gF.gf_irredPoly, 2));
      System.out.println("should get " + add.toString());
      System.out.println("But squared  is " + Poly.getRemainder(Poly.toPower(sqrt, gF.gf_irredPoly, 2), irredPoly, gF.gf_irredPoly).toString());
System.out.println("really sqrt? " + add.equals(Poly.getRemainder(Poly.toPower(sqrt, gF.gf_irredPoly, 2), irredPoly, gF.gf_irredPoly)));
      //calculate a and b s.t. (b * sqrt(x + synd)) = (a) mod (goppa poly)
      //where deg(a) <= floor(goppa deg / 2) and deg(b) <= floor((goppa deg - 1) / 2)
      Poly[] polys = Poly.partialGCDNew(sqrt, irredPoly, gF.gf_irredPoly);
      Poly polyA = polys[0];
      Poly polyB = polys[1];
//System.out.println("Poly A: " + polyA.toString());
//System.out.println("Poly B: " + polyB.toString());
      Poly test = Poly.multiply(polyB, sqrt, gF.gf_irredPoly);
//System.out.println("equal? " + Poly.getRemainder(test, irredPoly, gF.gf_irredPoly).toString() + " and " + polyA.toString());
      
      Poly errorLoc = Poly.multiply(polyA, polyA, gF.gf_irredPoly);
      polyB = Poly.multiply(polyB, polyB, gF.gf_irredPoly);
      polyB = Poly.multiply(Poly.xPoly, polyB, gF.gf_irredPoly);
     //errorLoc.coefficients ^= polyB.coefficients;
      errorLoc = Poly.add(errorLoc, polyB);
      Poly.getDegree(errorLoc);
System.out.println("Error Locator Polynomial: " + errorLoc.toString());
   
   
      
            
      Poly trace = McEliece.getTracePolynomial();

      List<gfPoly> roots = McEliece.runBerl(errorLoc, 0, trace);
      System.out.println("roots: ");
      for(int i = 0; i < roots.size(); i++)
      {
         System.out.println(roots.get(i).toString());
      }
      

      int[] decryptedError = new int[decoded.length];
      for(int i = 0; i < support.length; i++)
      {
         for(int k = 0; k < roots.size(); k++)
         {
            if(support[i].equals(roots.get(k)))
            {
               
               decryptedError[i] = 1;
            }
         }
      }


      int[][] errorMatrix = new int[decryptedError.length][1];
      
      
      for(int loc = 0; loc < decryptedError.length; loc++)
      {
         errorMatrix[loc][0] = decryptedError[loc];
      }
      //errorMatrix = multiply(invert(McEliece.permutation), errorMatrix);
      for(int eLoc = 0; eLoc < codeWord.length; eLoc++)
      {
         decryptedError[eLoc] = errorMatrix[eLoc][0];
      }
      System.out.println("error: ");
      McEliece.print(decryptedError);
      /*
      for(int i = 0; i < decoded.length; i++)
      {
         if(i < codeWord.length)
         {
            decoded[i] = codeWord[i] ^ decryptedError[i];
         }
         else
         {
            decoded[i] = decryptedError[i];
         }
      }
      */
      //codeWord = McEliece.add(codeWord, decryptedError);
System.out.println("After running decryption algorithm:");
McEliece.print(decoded);
      int[][] permInv = McEliece.invert(permutation);
      int[][] res = McEliece.multiply(permInv, vecTranspose(decoded));
      System.out.println("res: ");
      for(int i = 0; i < res.length; i++)
      {
         decoded[i] = res[i][0];
      }
      McEliece.print(decoded);
   
      System.out.println("Error perm? ");
      int[][] permErr = McEliece.multiply(permInv, errorMatrix);
      int[] permErrVec = new int[permErr.length];
      for(int i = 0; i < permErrVec.length; i++)
      {
         permErrVec[i] = permErr[i][0];
      }
            
      McEliece.print(permErrVec);
      int[] zeroes = new int[numErrors+1];
      int currGap = 0;
      int zeroesInd = 0;
      for(int errorInd = 0; errorInd < permErrVec.length; errorInd++)
      {
         if(permErrVec[errorInd] != 1)
         {
            currGap++;
         }
         else
         {
            zeroes[zeroesInd] = currGap;
            currGap = 0;
            zeroesInd++;
         }
      }
      zeroes[zeroesInd] = currGap;
      System.out.println("zeroes:");
      McEliece.print(zeroes);
      System.out.println("Decoded? ");
      McEliece.print(decodeCW(numErrors, permErrVec.length, zeroes));
      return decoded;
      
   }
   
   public static int[][] vecTranspose(int[] vec)
   {
      int[][] vecTransp = new int[vec.length][1];
      for(int loc = 0; loc < vec.length; loc++)
      {
         vecTransp[loc][0] = vec[loc];
      }
      
      return vecTransp;
   }
   
   public static Poly getSyndrome(int[] codeWord, int numErrors,
         int extDegree)
   {
      
      int currCoeffs = 0;
      int iter = 0; 
      System.out.println("in getSyn with ");
      McEliece.print(codeWord);
      Poly syndrome = new Poly(numErrors - 1);
      int[] synVec = codeWord;
      
      
      for(int index = synVec.length - 1; index >= 0; index--)
      {
         currCoeffs += (int)Math.pow(2, index % extDegree) * synVec[index];
         if(index % extDegree == 0)
         {
            syndrome.coefficients[iter] = new gfPoly(currCoeffs);
            System.out.println("gfPoly: " + syndrome.coefficients[iter].toString());
            currCoeffs = 0;
            iter++;
         }

      }
      
      /*
      
      for(int index = 0; index <= synVec.length-1; index++)
      {
         if(index % extDegree == 0 && index != 0)
         {
           
            syndrome.coefficients[iter] = new gfPoly(currCoeffs);
            currCoeffs = 0;
            iter++;
         }
         currCoeffs += (int)Math.pow(2, index % extDegree) * synVec[index];

        
        syndrome.coefficients[iter] = new gfPoly(currCoeffs);

      }
      */
      
      //coefficients go from highest to degree to lowest degree of x 
      //and lowest degree to highest degree of z
      return syndrome;
      
   }
   
   
}
