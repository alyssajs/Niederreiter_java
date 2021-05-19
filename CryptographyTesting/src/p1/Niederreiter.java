package p1;

import java.util.ArrayList;
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
      /*
      for(int test = 0; test < zeroes.length; test++)
      {
         System.out.print(zeroes[test] + " ");
      }
      System.out.println();
      */
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
      /*
      for(int t = 0; t < result.length; t++)
      {
         System.out.print(result[t] + " ");
      }
      System.out.println();
      */
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
   
   public static int[] encrypt(int[] message, int extDegreeIn, int numErrorsIn, int suppSize)
   {
	   do
	   {
	      parityCheck = getParityCheck(extDegreeIn, numErrorsIn, suppSize);
	   }while(!isSystematic(parityCheck, parityCheck[0].length-parityCheck.length));
   
	  // System.out.println("Extension Degree: " +  extDegree);
	   
	  // System.out.println("Number of Errors: " + numErrors);
	   
	   /*
	   System.out.println("support: ");
	   print(support);
	   */
	   
	   System.out.println("Goppa polynomial: " + irredPoly.toString());
	
	   
	  System.out.println("Field polynomial: " + gF.gf_irredPoly);
	   
	  System.out.println("Parity Check: ");
	 print(parityCheck);
	  
	   
	   	      
	   
	   
	   
      //placeholder matrices to calculate invertible matrix
      int dim = parityCheck.length;
      int[][] mA = new int[dim][dim];
      int[][] mT = new int[dim][dim];
      
      //calculate random invertible matrix 
      scramble = getRandMatrix(mA, mT, dim);
     
      
      //calculate permutation matrix 
      permutation = createPermMatrix(parityCheck[0].length);
     
      //dELETE this 
      int[] errorVec = new int[parityCheck[0].length];
      gfPoly[] pOne = new gfPoly[2];
      gfPoly[] pTwo = new gfPoly[2];
      int numOnes = 0;
      for(int i = 0; i < message.length; i++)
      {
    	  if(message[i] == 1)
    	  {
    		  errorVec[i] = 1;
    	  }
      }
      int[][] errorVecMat = multiply(permutation, vecTranspose(errorVec));
      System.out.println("permuted error: ");
      print(errorVecMat);
      System.out.println("permInv? ");
      print(multiply(invert(permutation), errorVecMat));
      System.out.println("errorVecMat dims: " + errorVecMat.length + " x " + errorVecMat[0].length);
      for(int i = 0; i < errorVecMat.length; i++)
      {
    	  if(errorVecMat[i][0] == 1)
    	  {
    		  if(numOnes == 0)
    		  {
    			  pOne[0] = support[i];
    			  numOnes++;
    		  }
    		  else
    		  {
    			  //pTwo[0] = support[i];
    		  }
    	  }
    		  
      }
      
      pOne[1] = new gfPoly(1);
      //pTwo[1] = new gfPoly(1);
      Poly polyOne = new Poly(pOne);
      /*
      Poly polyTwo = new Poly(pTwo);
      System.out.println("polyOne coeff vector:");
      print(pOne);
      System.out.println("polyTwo coeff vector: ");
     print(pTwo);
      System.out.println("polyOne coeffs 0 null? " + pOne[0]==null);
      System.out.println("polyTwo coeffs null? " + pTwo[0]==null);
     // System.out.println("polyTwo.coeffs null? " + polyOne.coefficients==null);
      System.out.println("polyTwo.coeffs null? " + polyTwo.coefficients==null);

      System.out.println("Two polys are " + polyOne.toString() + " and " + polyTwo.toString());
      */
      System.out.println("errorLoc should be: " + polyOne.toString());
      /*
      Poly prod = Poly.multiply(polyOne, polyTwo, gF.gf_irredPoly);
      System.out.println("error loc should be: " + prod.toString());
      Poly trace = getTracePolynomial();
      List<gfPoly> theseRoots = runBerl(prod, 0, trace);
      System.out.println("these roots: ");
      for(int i = 0; i < theseRoots.size(); i++)
      {
    	  System.out.println(theseRoots.get(i) + ", ");
      }
      System.out.println("^^^^^^^");
      */
      //^^^^^^^^^^^^^^^^^
      int[][] messageMatTransp = new int[message.length][1];
      
      
      for(int matRow = 0; matRow < message.length; matRow++)
      {
         messageMatTransp[matRow][0] = message[matRow];
      }
      
      //multiply S * parityCheck * permutation * messageTranspose
      int[][] encryptedMat = multiply(scramble, parityCheck);

      encryptedMat = multiply(encryptedMat, permutation);
System.out.println("Public key:");
print(encryptedMat);
      encryptedMat = multiply(encryptedMat, messageMatTransp);
      
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
   
   /*
   static int[][] createParityCheck(int extDegreeIn, int numErrorsIn, int suppSize)
   {
	   int[][] parityCheckMatrix;
	   do
	   {
	      parityCheckMatrix = getParityCheck(extDegreeIn, numErrorsIn, suppSize);
	   }while(!isSystematic(parityCheckMatrix, parityCheckMatrix[0].length-parityCheckMatrix.length));
	   parityCheck = parityCheckMatrix;

      extDegree = extDegreeIn;
      numErrors = numErrorsIn;
      support = McEliece.support;
      gF = McEliece.gF;
 System.out.println("fieldPoly: " + gF.gf_irredPoly);
 System.out.println("irredPoly: " + irredPoly.toString());
      polyParity = McEliece.polyParity;  
      unreducedParity = McEliece.unreducedParity;
      return parityCheck;
   }
   */
   static int[] decrypt(int[] ciphertext, int[] actual)
   {
      //multiply by random matrix inverse
      int[][] randInverse = invert(scramble);
     
      /*
      int[][] decryptedTemp = new int[1][ciphertext.length];
      
      
      decryptedTemp[0] = ciphertext;
      
      decryptedTemp = multiply(decryptedTemp, randInverse);
      */
       
      int[][] decryptedTemp = vecTranspose(ciphertext);
      decryptedTemp = multiply(randInverse, decryptedTemp);
      
      System.out.println("after multiplying:");
     print(decryptedTemp);
    
      
      int[][] decryptedNext = new int[1][parityCheck[0].length];
    
      //use syndrome decoder (patterson/BTA??)
      
      int[] unscrambled = new int[decryptedTemp.length];
      for(int loc = 0; loc < unscrambled.length; loc++)
      {
         unscrambled[loc] = decryptedTemp[loc][0];
      }
      System.out.println("passing in for decryption:");
      print(unscrambled);
      decryptedNext[0] = pattersonDecode(unscrambled, parityCheck, actual);
      
      /*
      int[][] decryptedNext = new int[1][parityCheck[0].length];
      decryptedNext[0] = pattersonDecode(ciphertext, parityCheck);
      */
      //then multiply by permInverse
      int[][] permInverse = invert(permutation);
      
      decryptedNext = multiply(decryptedNext, permInverse);
      
      return decryptedNext[0];
   }
   

   public static int[] pattersonDecode(int[] codeWord, int[][] parityCheck, int[] actual)
   {
      System.out.println("codeWord.length+ " + codeWord.length);
      System.out.println("support size " + support.length);
      
      System.out.println("in here with: ");
      print(codeWord);
      //int[][] scrambleInv = invert(scramble);

  
      
     // int[][] codeTemp = new int[1][codeWord.length];
      //codeTemp[0] = codeWord;
      //System.out.println("codeTemp dims " + codeTemp.length + ", " + codeTemp[0].length);
     // System.out.println("sI dims: " + scrambleInv.length + ", " + scrambleInv[0].length);
      //codeTemp = multiply(codeTemp, scrambleInv);

      //codeWord = codeTemp[0];
      
      int[] decoded = new int[parityCheck[0].length];
   
      //Multiply ciphertext by parity check matrix to get syndrome polynomial
      //Poly syndrome = getSyndrome(codeWord, numErrors, extDegree);
      Boolean gotAnswer = false;

      Poly syndrome = altSyndrome(codeWord);
      System.out.println("Syndrome: " + syndrome.toString());  
      //gfPoly inv = gfPoly.getModularInverse(syndrome.coefficients[1], gF.gf_irredPoly);
      //syndrome = Poly.multiply(syndrome, inv, gF.gf_irredPoly);
      //System.out.println("now syndrome is: " + syndrome.toString());
      /*
      gfPoly pOne = new gfPoly(7);
      gfPoly pTwo = new gfPoly(1);
      gfPoly pThree = new  gfPoly(1);
      syndrome = new Poly(new gfPoly[] {pOne, pTwo, pThree});
      */
      while(!gotAnswer)
      {
      syndrome = Poly.getRandSyn(numErrors-1, extDegree);
      System.out.println("randomly generated syndrome: " + syndrome.toString());
  

      //Invert syndrome polynomial mod goppa polynomial 
      Poly synInverse = Poly.getModularInverse(syndrome, irredPoly, gF.gf_irredPoly).get(1);
System.out.println("really inv? " + Poly.multiply(synInverse, syndrome, gF.gf_irredPoly, irredPoly).toString());
      //add f(x) = x to inverse
      //Polynomial add = new Polynomial(synInverse.coefficients + 2);
      Poly add = Poly.add(synInverse, Poly.xPoly);
      //calculate sqrt(x + syndrome inverse)
//System.out.println("trying to get sqrt of " + add.toString());
//System.out.println("polys are " + irredPoly.toString() + " and " + gF.gf_irredPoly);
      Poly sqrt = Poly.calcSqrt(add, irredPoly, gF.gf_irredPoly.degree, gF.gf_irredPoly);
//      System.out.println("got sqrt: " + sqrt.toString());
// System.out.println("sqrt is: " + sqrt.toString());
      //System.out.println("before reducing: " + Poly.toPower(sqrt, gF.gf_irredPoly, 2));
//      System.out.println("should get " + add.toString());
//      System.out.println("But squared  is " + Poly.getRemainder(Poly.toPower(sqrt, gF.gf_irredPoly, 2), irredPoly, gF.gf_irredPoly).toString());
//System.out.println("really sqrt? " + add.equals(Poly.getRemainder(Poly.toPower(sqrt, gF.gf_irredPoly, 2), irredPoly, gF.gf_irredPoly)));
      //calculate a and b s.t. (b * sqrt(x + synd)) = (a) mod (goppa poly)
      //where deg(a) <= floor(goppa deg / 2) and deg(b) <= floor((goppa deg - 1) / 2)
      Poly[] polys = Poly.partialGCDNew(sqrt, irredPoly, gF.gf_irredPoly);
      Poly polyA = polys[0];
      Poly polyB = polys[1];
//System.out.println("Poly A: " + polyA.toString() + " should have degree smaller than " + (numErrors)/2);
//System.out.println("Poly B: " + polyB.toString() + " should have degree smaller than " + (numErrors-1)/2);
      Poly test = Poly.multiply(polyB, sqrt, gF.gf_irredPoly);
//System.out.println("equal? " + Poly.getRemainder(test, irredPoly, gF.gf_irredPoly).toString() + " and " + polyA.toString());
      
      Poly errorLoc = Poly.multiply(polyA, polyA, gF.gf_irredPoly);
      polyB = Poly.multiply(polyB, polyB, gF.gf_irredPoly);
      polyB = Poly.multiply(Poly.xPoly, polyB, gF.gf_irredPoly);
      errorLoc = Poly.add(errorLoc, polyB);
      Poly.getDegree(errorLoc);
System.out.println("Error Locator Polynomial: " + errorLoc.toString());
//gfPoly inv = gfPoly.getModularInverse(errorLoc.coefficients[2], gF.gf_irredPoly);
//System.out.println("monic errorLoc " + Poly.multiply(errorLoc, inv, gF.gf_irredPoly).toString());   
   
      
            
      Poly trace = getTracePolynomial();

      List<gfPoly> roots = runBerl(errorLoc, 0, trace);
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
      print(decryptedError);
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
      //codeWord = add(codeWord, decryptedError);
System.out.println("After running decryption algorithm:");
print(decoded);
      int[][] permInv = invert(permutation);
      int[][] res = multiply(permInv, vecTranspose(decoded));
      System.out.println("res: ");
      for(int i = 0; i < res.length; i++)
      {
         decoded[i] = res[i][0];
      }
      print(decoded);
   
      System.out.println("Error perm? ");
      int[][] permErr = multiply(permInv, errorMatrix);
      int[] permErrVec = new int[permErr.length];
      for(int i = 0; i < permErrVec.length; i++)
      {
         permErrVec[i] = permErr[i][0];
      }
            
      print(permErrVec);
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
      print(zeroes);
      System.out.println("Decoded? ");
      print(decodeCW(numErrors, permErrVec.length, zeroes));
      
      if(permErrVec[1] == 1)
      {
    	  gotAnswer = true;
      }
      }
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
   public static gfPoly[][] vecTranspose(gfPoly[] vec)
   {
      gfPoly[][] vecTransp = new gfPoly[vec.length][1];
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
      int iter = numErrors-1; 
      System.out.println("in getSyn with ");
      print(codeWord);
      Poly syndrome = new Poly(numErrors - 1);
      int[] synVec = codeWord;
      gfPoly[] synPolys = new gfPoly[codeWord.length/extDegree];
      
      
      
      for(int index = 0; index <= synVec.length-1; index++)
      {
         if(index % extDegree == 0 && index != 0)
         {
           
            syndrome.coefficients[iter] = new gfPoly(currCoeffs);
            synPolys[iter] = new gfPoly(currCoeffs);

            currCoeffs = 0;
            iter--;
         }
         currCoeffs += (int)Math.pow(2, index % extDegree) * synVec[index];

        
        syndrome.coefficients[iter] = new gfPoly(currCoeffs);

      }
      /*
      
      for(int index = synVec.length - 1; index >= 0; index--)
      {
         currCoeffs += (int)Math.pow(2, index % extDegree) * synVec[index];
         if(index % extDegree == 0 && index != 0)
         {
            syndrome.coefficients[iter] = new gfPoly(currCoeffs);
            synPolys[iter] = new gfPoly(currCoeffs);
            System.out.println("gfPoly: " + syndrome.coefficients[iter].toString());
            currCoeffs = 0;
            iter++;
            
            
         }
         
      }
      */
      
      syndrome.coefficients[iter] = new gfPoly(currCoeffs);
      
      /*
      synPolys[iter] = new gfPoly(currCoeffs);

      gfPoly[][] coeffMat = new gfPoly[numErrors][numErrors];
      for(int row = 0; row < numErrors; row++)
      {
    	  for(int col = 0; col < numErrors; col++)
    	  {
    		  if(col >= row) 
    		  {
    			  coeffMat[row][col] = irredPoly.coefficients[numErrors - (col - row)];
    		  }
    		  else
    		  {
    			  coeffMat[row][col] = new gfPoly(0);
    		  }
    	  }
      }
      McEliece.print(coeffMat);
      McEliece.print(synPolys);
      gfPoly[][] prod = McEliece.multiply(coeffMat, vecTranspose(synPolys), gF.gf_irredPoly);
      McEliece.print(prod);
      gfPoly[] prodVec = new gfPoly[prod.length];
      for(int i = 0; i < prodVec.length; i++)
      {
    	  prodVec[i] = prod[i][0];
      }
      syndrome = new Poly(prodVec);
      */
      


      

      //coefficients go from highest to degree to lowest degree of x 
      //and lowest degree to highest degree of z
      return syndrome;
      
   }
   public static void print(int[] vec)
	{
	   System.out.print("[");
	   
	   for(int i = 0; i < vec.length; i++)
	   {
	      System.out.print(vec[i]);
	      if(i != vec.length - 1)
	      {
	         System.out.print(", ");
	      }
	      
	     
	   }
	   
	   System.out.print("] \n");
	}
	
	public static void print(gfPoly[] vec)
	{
	   System.out.print("[");
     
     for(int i = 0; i < vec.length; i++)
     {
        System.out.print(vec[i].toString());
        if(i != vec.length - 1)
        {
           System.out.print(", ");
        }
         
       
     }
     
     System.out.print("] \n");
	}

	
	/**
	 * Creates a random length x length permutation matrix
	 * @param length Length of input messages
	 * @return Random permutation matrix of dimensions length x length
	 */
	public static int[][] createPermMatrix(int length)
	{
		int[][] permMatrix = new int[length][length];
		
		//identity matrix rows that permMatrix rows will be swapped with
		int[] swappingWith = new int[length];
		
		int adding;
		
		for(adding = 0; adding < length; adding++)
		{
			swappingWith[adding] = length;
		}
		
		int[][] idMatrix = createIDMatrix(length);
		
		int toSwap, index, swapping;
		
		for(index = 0; index < length; index++)
		{
			do
			{
				toSwap = (int)(Math.random() * length);
			
			}while(isInArray(swappingWith, toSwap));
			
			swappingWith[index] = toSwap;
		}
		
		for(swapping = 0; swapping < length; swapping++)
		{
			permMatrix[swapping] = idMatrix[swappingWith[swapping]];
			
		}
		
		return permMatrix;
	}
	
	public static boolean isInArray(int[] array, int testVal)
	{
		for(int i = 0; i < array.length; i++)
		{
			if(array[i] == testVal)
			{
				return true;
			}
		}
		
		return false;
	}
	
	public static int[][] createIDMatrix(int length)
	{
		int index;
		int[][] idMatrix = new int[length][length];
		
		for(index = 0; index < length; index++)
		{
		
			idMatrix[index][index] = 1;
		
		}
		
		return idMatrix;
	}
	

	public static void print(int[][] array)
	{
		int rowIndex, columnIndex;
		
		for(rowIndex = 0; rowIndex < array.length; rowIndex++)
		{
		   System.out.print("[");
			for(columnIndex = 0; columnIndex < array[0].length; columnIndex++)
			{
			   
			   System.out.print(array[rowIndex][columnIndex]);
			   if(columnIndex != (array[0].length - 1))
			   {
			      System.out.print(", ");
			            
			   }
			   
			}
			System.out.print("]");
			
			System.out.println();
		}
	}
	  
   public static void print(Poly[][] array)
   {
      int rowIndex, columnIndex;
      
      for(rowIndex = 0; rowIndex < array.length; rowIndex++)
      {
         for(columnIndex = 0; columnIndex < array[0].length; columnIndex++)
         {
            System.out.print(array[rowIndex][columnIndex].toString() + " " );
         }
         
         System.out.println();
      }
   }
	
	public static int[][] multiply(int[][] matrix1, int[][] matrix2)
	{
		if(matrix1[0].length != matrix2.length)
		{
			return null; 
		}
		int[][] result = new int[matrix1.length][matrix2[0].length];
				
		
		
		int matrix1Row, matrix2Column, iterIndex;
		
		for(matrix1Row = 0; matrix1Row < matrix1.length; matrix1Row++)
		{
			for(matrix2Column = 0; matrix2Column < matrix2[0].length;
					matrix2Column++)
			{
				result[matrix1Row][matrix2Column] = 0;
				for(iterIndex = 0; iterIndex < matrix2.length;
						iterIndex++)
					{
						result[matrix1Row][matrix2Column] += 
								(matrix1[matrix1Row][iterIndex] *
								matrix2[iterIndex][matrix2Column]) ;

					}
				result[matrix1Row][matrix2Column] = result[matrix1Row][matrix2Column] % 2;
				
			}
				
		}
				
		return result;
	}
	
	public static int[][] getRandMatrix(int[][] mA, int[][] mT, int dim)
	{
	   ArrayList<int[][]> matrices = genMatrix(mA, mT, dim);
	   return multiply(matrices.get(0), matrices.get(1));
	}

	public static ArrayList<int[][]> genMatrix(int[][] mA, int[][] mT, int dim)
	{
		
		
		if(dim == 1)
		{
			mA[0][0] = 1;
			
			mT[0][0] = 1;
			
		}
		else if(dim > 1)
		{
			int[] vector = new int[dim];
		
			int vecIndex, firstNonzero;
		
			do
			{
				for(vecIndex = 0; vecIndex < vector.length; vecIndex++)
				{
					vector[vecIndex] = (int)(Math.random() * 2);
				}
				
			}	
			while(isZeroVector(vector));
		
			vecIndex = 0;
		
			firstNonzero = vector.length;
		
			while(firstNonzero == vector.length)
			{
				if(vector[vecIndex] != 0)
				{
					firstNonzero = vecIndex;
				}
			
				vecIndex++;
			}
			
			
			int[] firstRow = new int[dim];
			firstRow[firstNonzero] = 1;
			
			mA[0] = firstRow;
			
			
			mT[firstNonzero] = vector;
			
		
			int copyingRow, copyingColumn;
			
			int newRow = 0;
			
			int newColumn = 0;
         //create the minors

			int[][] newMatrixA = new int[mA.length - 1][mA[0].length - 1];
						
			int[][] newMatrixT = new int[mT.length - 1][mT[0].length  - 1];
			
			newRow = 0;
			newColumn = 0;
			copyingRow = 0;
			copyingColumn = 0;
		

			ArrayList<int[][]> returned = genMatrix(newMatrixA, newMatrixT, dim-1); 
			newMatrixA = returned.get(0);
			for(newRow = 1; newRow < mA.length; newRow++)
			{
			      for(newColumn = 0, copyingColumn = 0; newColumn < mA[0].length; newColumn++)
			      {
			         if(newColumn != firstNonzero)
			         {
			      
			            mA[newRow][newColumn]= newMatrixA[copyingRow][copyingColumn];

			            copyingColumn++;
			         }

			      }
			      copyingRow++;
			}
			
			

			newMatrixT = returned.get(1);
	      copyingRow = copyingColumn = 0;
         for(newRow = 0; newRow < mT.length; newRow++)
         {
            
            if(newRow != firstNonzero)
            {
               for(newColumn = 0, copyingColumn = 0; newColumn < mT[0].length; newColumn++)
               {
                  if(newColumn != firstNonzero)
                  {

                     mT[newRow][newColumn]= newMatrixT[copyingRow][copyingColumn];
                     copyingColumn++;
                  }
               }
               copyingRow++;  
            }
         }



		}
		
		ArrayList<int[][]> arrays = new ArrayList<int[][]>();
		
		arrays.add(mA);
		
		arrays.add(mT);
		return arrays;
		
	}
	

	
	public static boolean isZeroVector(int[] vector)
	{
		int index;
		
		for(index = 0; index < vector.length; index++)
		{
			if(vector[index] != 0)
			{
				return false;
			}
		}
		
		return true;
	}
	
	public static int[][] invert(int[][] matrix)
	{
//System.out.println("inverting:");
//print(matrix);
	   int[][] orig = matrix;
	   int[][] newMatrix = new int[matrix.length][2 * matrix.length];
	   for(int addOne = 0; addOne < matrix.length; addOne++)
	   {
	      newMatrix[addOne][matrix.length + addOne] = 1;
	   }
  
	 
	   int copyingRow, copyingColumn;
	   
	   for(copyingRow = 0; copyingRow < matrix.length; copyingRow++)
	   {
	      for(copyingColumn = 0; copyingColumn < matrix.length; copyingColumn++)
	      {
	         newMatrix[copyingRow][copyingColumn] = matrix[copyingRow][copyingColumn];
	      }
	   }
	   
	   matrix = newMatrix;
	   //System.out.println("row reducing:");
	  //print(matrix);  
	   
	   int searchRow;
	   
	  for(int j = 0; j < matrix.length; j++)
	  {
	     
	      searchRow = j;
	      
	      boolean found = false;
	      /*
	      while(found == false && k < matrix.length)
	      {
	         if(matrix[k][j] == 1)
	           {
	             searchRow = k;
	             found = true;
	           } 
	         k++;
	      }
	      */
        
	     for(int k = searchRow; k < matrix.length; k++)
	     {
	        if(matrix[k][j] == 1)
	        {
	          searchRow = k;
	        }
	       
	     }
	     
	      
	     if(matrix[searchRow][j] == 1)
        {
           int[] temp = matrix[searchRow];
           matrix[searchRow] = matrix[j];
           matrix[j] = temp;
            /*          
           for(int u = j + 1; u < matrix.length; u++)
           {
              for(int matrCol = 0; matrCol < matrix[0].length; matrCol++)
              {
System.out.println("multiplying " + matrix[u][j] + " which is at " + u + ", " + j + " by entry at " + j + ", " + matrCol + " and assigning to " + u + ", " + matrCol);                
                 matrix[u][matrCol] ^= (matrix[u][j] * matrix[j][matrCol]);
              }
           }
           System.out.println("now matrix is: ");
           print(matrix);
            */
        }

	     
	       int clearingCol = j;
	        for(int clearingRow = clearingCol + 1; clearingRow < matrix[0].length/2; clearingRow++)
	        {
	           if(matrix[clearingRow][clearingCol] == 1)
	           {
	              matrix[clearingRow] = add(matrix[clearingRow], matrix[clearingCol]);
	           }
	        }

	     
	
	  }
	  //System.out.println("after clearing below");
	 // print(matrix);
	  int pivotCol = matrix[0].length / 2 - 1;
	  int pivotRow = matrix.length - 1;
	  
	  int currentRow, currentCol;
	  
     for(currentCol = pivotCol; currentCol >= 0; currentCol--)
	  {
        for(currentRow = matrix.length - 2; currentRow >= 0; currentRow--)
	     {
	        if(matrix[currentRow][currentCol] == 1 && (currentRow != pivotRow))
	        {
	           for(int addCol = 0; addCol < matrix[0].length; addCol++)
	           {
	              matrix[currentRow][addCol] ^= matrix[pivotRow][addCol];
	           }
	           
	           
	        }
           pivotCol--;

	     }
	     
        pivotRow--;

	  }
	  

      int[][] matrixInverse = new int[matrix.length][matrix[0].length / 2];
      for(int inverseRow = 0; inverseRow < matrix.length; inverseRow++)
      { 
         for(int inverseColumn = matrix.length; inverseColumn < 2 *matrix.length
               ;inverseColumn++)
         {
            matrixInverse[inverseRow][inverseColumn - matrix.length]
                  = matrix[inverseRow][inverseColumn];
         }  
      }
      
      
     
      matrix = matrixInverse;

      return matrix;
      
      
      
	}
	
	public static int[][] getParityCheck(int extenDegree, int numErrorsIn, int suppSize)
	{
       gF = new GaloisField(extenDegree);
       extDegree = extenDegree;
      numErrors = numErrorsIn;
      irredPoly = Poly.getIrredPoly(numErrors, gF.gf_irredPoly);
      support = gF.getSupport(irredPoly, suppSize);
      /*
      System.out.println("support:");
      print(support);
      */

      
	   Poly[][] mat1 = new Poly[irredPoly.degree][irredPoly.degree];
	   Poly[][] mat2 = new Poly[irredPoly.degree][support.length];
	   Poly[][] mat3 = new Poly[support.length][support.length];

	   
	   //matrix 1 is a lower triangular matrix with the coefficients of the irredPoly
	   for(int i = 0; i < mat1.length; i++)
	   {
	      for(int j = 0; j <= i; j++)
	      {
	         if(j == i)
	         { 
	            //mat1[i][j] = new Polynomial(1);
	            mat1[i][j] = Poly.onePoly;
	         }
	         else
	         {
	           // mat1[i][j] = irredPoly.coefficients[irredPoly.degree - i];
	            //TODO: correct change?
	           //mat1[i][j] = new Polynomial(irredPoly.coefficients >> (irredPoly.degree - (i -j)) & 1);
	            mat1[i][j] = new Poly(new gfPoly[] {irredPoly.coefficients[irredPoly.degree - (i-j)]}) ;
	         }

	      }
	      for(int j = i + 1; j < mat1[0].length; j++)
         {
            //mat1[i][j] = new Polynomial(0);
	         mat1[i][j] = Poly.zeroPoly;

         }
	      

	   }

	   //matrix2 has the powers of each element of the support
	   for(int i = 0; i < mat2.length; i++)
	   {
	      for(int j = 0; j < mat2[0].length; j++)
	      {
	         mat2[i][j] = new Poly(gfPoly.toPower(support[j], gF.gf_irredPoly, i));
	         
	      }
	   }

	   //matrix3 has inverses of g(support) on diagonal 
	   for(int i = 0; i < mat3.length; i++)
	   {
	      for(int j = 0; j < mat3[0].length; j++)
	      {
	         if(i == j)
	         {  
	            //TODO: mod out by # of errors polynomial or galois field polynomial?
	            gfPoly entry = Poly.eval(irredPoly, support[i], gF.gf_irredPoly);
	            
	            
	           //entry = Polynomial.getRemainder(entry, gF.gf_irredPoly);
	            
	            mat3[i][j] = new Poly(gfPoly.getModularInverse(entry, gF.gf_irredPoly));

	    

	         }
	         else
	         {
	            //mat3[i][j] = new Polynomial(0);
	            mat3[i][j] = Poly.zeroPoly;
	         }
	      }
	   }

	   
      Poly[][] polyMatr = multiply(multiply(mat1, mat2), mat3);
      polyParity = polyMatr;

	   
	   int pcRow = 0;
	   int pcCol = 0;
	   

	   
	   int[][] parityCheckMatrix = new int[irredPoly.degree * extDegree][support.length];
	   for(int i = 0; i < polyMatr.length; i++)
	   {
	      for(int j = 0; j < polyMatr[0].length; j++)
	      {
	         //pcRow = irredPoly.degree*i?
	         pcRow = gF.gf_irredPoly.degree*i;
	         pcCol = j;
	         for(int k = 0; k < gF.gf_irredPoly.degree; k++)
	         {
	            if(k <= polyMatr[i][j].coefficients[0].degree)
	            {
	               parityCheckMatrix[pcRow][pcCol] = (polyMatr[i][j].coefficients[0].coeffs >> k) & 1;
	            }
	            else
	            {
	               parityCheckMatrix[pcRow][pcCol] = 0;
	            }
	            pcRow++;
	         }
	      }
	   }

      unreducedParity=new int[parityCheckMatrix.length][parityCheckMatrix[0].length];
      for(int row = 0; row < unreducedParity.length; row++)
      {
         for(int col = 0; col < unreducedParity[0].length; col++)
         {
            unreducedParity[row][col] = parityCheckMatrix[row][col];
         }
      }
   
	   parityCheckMatrix = rowReduce(parityCheckMatrix, parityCheckMatrix[0].length - parityCheckMatrix.length);

	   return parityCheckMatrix;
	}
	
	public static Poly[][] multiply(Poly[][] matrix1, Poly[][] matrix2)
   {
      if(matrix1[0].length != matrix2.length)
      {
         return null; 
      }
      Poly[][] result = new Poly[matrix1.length][matrix2[0].length];
            
      
      
      int matrix1Row, matrix2Column, iterIndex;
      
      for(matrix1Row = 0; matrix1Row < matrix1.length; matrix1Row++)
      {
         for(matrix2Column = 0; matrix2Column < matrix2[0].length;
               matrix2Column++)
         {
            //result[matrix1Row][matrix2Column] = new Polynomial(0);
            result[matrix1Row][matrix2Column] = Poly.zeroPoly;
            for(iterIndex = 0; iterIndex < matrix2.length;
                  iterIndex++)
               {
                  result[matrix1Row][matrix2Column] = Poly.add(result[matrix1Row][matrix2Column], 
                        Poly.multiply(matrix1[matrix1Row][iterIndex],
                        matrix2[iterIndex][matrix2Column], gF.gf_irredPoly)) ;

               }

            
         }
            
      }
            
      return result;
   }
	
	public static Poly[][] reduce(Poly[][] matr, Poly irredPoly)
	{
	   for(int i = 0; i < matr.length; i++)
	   {
	      for(int j = 0; j < matr[0].length; j++)
	      {
	        matr[i][j] = Poly.getRemainder(matr[i][j], irredPoly, gF.gf_irredPoly);
	      }
	   }
	   
	   return matr;
	}
	
	public static int[][] rowReduce(int[][] matrix, int startingCol)
	{
	   int searchRow;
	   
	     for(int j = startingCol; j < matrix[0].length; j++)
	     {        
	        //find a 1 in column j
	        searchRow = j - startingCol; //???
	        for(int k = j - startingCol; k < matrix.length; k++)
	        {
	      
	           if(matrix[k][j] == 1)
	           {
	             searchRow = k;
	           }
	          
	        }

	        //if we found one, swap it so it's in the right position 
	        if(matrix[searchRow][j] == 1)
	        {

	           int[] temp = matrix[searchRow];
	           matrix[searchRow] = matrix[j - startingCol];
	           matrix[j - startingCol] = temp;
	            //zero everything below it?          
	             for(int i = j - startingCol +1; i < matrix.length; i++)
	             {
	                if(matrix[i][j] == 1)
	                {
	                   for(int add = 0; add < matrix[0].length; add++)
	                   {
	                      matrix[i][add] ^= matrix[j - startingCol][add];
	                   }
	                }
	             }
	            
	        }
	     }


	    
	     /* 
	     for(currentCol = pivotCol; currentCol >= 0; currentCol--)
	     {
	        for(currentRow = matrix.length - 2; currentRow >= 0; currentRow--)
	        {
	           if(matrix[currentRow][currentCol] == 1 && (currentRow != pivotRow))
	           {
	              for(int addCol = 0; addCol < matrix[0].length; addCol++)
	              {
	                 matrix[currentRow][addCol] ^= matrix[pivotRow][addCol];
	              }
	              
	              
	           }
	           pivotCol--;

	        }
	        
	        pivotRow--;

	     }
	     */

	     //reduced row echelon form
	     
	     for(int matrCol = startingCol; matrCol < matrix[0].length; matrCol++)
	     {
	        for(int matrRow = 0; matrRow < matrCol - startingCol; matrRow++)
	        {
	           //if upper entries are 1's, add the proper row to make them zero
	           if( matrix[matrRow][matrCol] == 1)
	           {
	              //add row matrCol to matrRow??
	              for(int i = 0; i < matrix[0].length; i++)
	              {
	                 matrix[matrRow][i] ^= matrix[matrCol - startingCol][i];
	              }
	           }
	        }
	     }
	     
	     
	     return matrix;
	     
	}
	
	
	public static int[][] getTranspose(int[][] matrix)
	{
	   int[][] transpose = new int[matrix[0].length][matrix.length];
	   
	   for(int col = 0; col < matrix[0].length; col++)
	   {
	      for(int row = 0; row < matrix.length; row++)
	      {
	         transpose[col][row] = matrix[row][col];
	      }
	   }
	   
	   return transpose;
	}
	

	static int[] add(int[] vec1, int[] vec2) 
	{
      if(vec1.length == vec2.length)
      {
         int[] result = new int[vec1.length];
         for(int loc = 0; loc < vec1.length; loc++)
         {
            result[loc] = vec1[loc] ^ vec2[loc];
         }
         return result;
      }
      return null;
   }

	   public static List<gfPoly> runBerl(Poly errorLoc, int basisLoc, Poly trace)
		{

		   if(errorLoc.degree <= 1)
		   {
		      List<gfPoly> roots = Poly.getRoot(errorLoc, gF.gf_irredPoly);
		      return roots;
		   }
		   gfPoly basis_i = gfPoly.toPower(new gfPoly("10", 1), gF.gf_irredPoly,
		         basisLoc);
		   
		   Poly evaluated = Poly.eval(trace, Poly.multiply(Poly.xPoly, basis_i, 
		         gF.gf_irredPoly), gF.gf_irredPoly);
		   

		   Poly sig_0 = Poly.gcd(errorLoc, evaluated, gF.gf_irredPoly);
		   
		   Poly eval2 = Poly.add(evaluated, Poly.onePoly);
		   
		   
		   Poly sig_1 = Poly.gcd(errorLoc, eval2, gF.gf_irredPoly);
		
		   List<gfPoly> result1 = runBerl(sig_0, basisLoc + 1, trace);
		   List<gfPoly> result2 = runBerl(sig_1, basisLoc + 1, trace);
		   for(int i = 0; i < result2.size(); i++)
		   {
		      result1.add(result2.get(i));
		   }
		   return result1;
		   
		}
	   
		public static Poly getTracePolynomial()
		{
		   Poly trace = new Poly((int)Math.pow(2, extDegree-1));
		   for(int i = 1; i <= trace.degree; i *= 2)
		   {
		      trace.coefficients[i] = gfPoly.gfOne;
		   }
	      
		   return trace;
		}
		public static boolean isSystematic(int[][] matr, int startingCol)
		{
		   for(int i = 0; i < matr.length; i++)
		   {
		      for(int j = startingCol; j < matr[0].length; j++)
		      {
		         if(i == j - startingCol)
		         {
		            if(matr[i][j] != 1)
		            {
		               return false;
		            }
		         }
		      }
		   }
		   
		   return true;
		}
		

		public static Poly altSyndrome(int[] encoded)
		{
			/*
			System.out.println("in here with irredPoly: " + irredPoly.toString());
			System.out.println("and gfPoly is " + gF.gf_irredPoly);
			System.out.println("and cw is ");
			McEliece.print(encoded);
			McEliece.print(support);
			*/
		   Poly factor;
		   
		   Poly syndrome = new Poly(Poly.zeroPoly, 0);
		   for(int i = 0; i < encoded.length; i++)
		   {
		      if(encoded[i] == 1)
		      {
		    	  System.out.println("Getting inverse of x - " + support[i].toString());
		         factor = Poly.getModularInverse(Poly.add(Poly.xPoly, new Poly(support[i])), irredPoly, gF.gf_irredPoly).get(1);
		         syndrome = Poly.add(syndrome,  factor);
		         System.out.println("now syndrome is: " + syndrome.toString());
		      }
		      
		      
		   }
		   
			/*
			System.out.println("polyParity: ");
			print(polyParity);
			
			System.out.println("extending: ");
			print(encoded);
			int[] encCopy = new int[polyParity[0].length];
			for(int i = 0; i < encoded.length; i++)
			{
				encCopy[i] = encoded[i];
			}
			System.out.println("extended: ");
			print(encCopy);
			System.out.println("mult");
			gfPoly[][] res = (mult(encCopy, polyParity, gF.gf_irredPoly));
		
			   Poly syndrome = new Poly();
			   Poly newFact;
		      //if code word vector is 1, then add the factor	  
		      for(int parityCol = 0; parityCol < polyParity[0].length; parityCol++)
		      {
		         if(encCopy[parityCol] == 1)
		         {
		            for(int parityRow = 0; parityRow < polyParity.length; parityRow++)
		            {

		  
		               newFact = Poly.toPower(Poly.xPoly, gF.gf_irredPoly, polyParity.length - parityRow - 1);
		               newFact = Poly.multiply(newFact, polyParity[parityRow][parityCol], gF.gf_irredPoly);
		               syndrome = syndrome.add(syndrome, newFact);
		            }
		         } 
		         
		      }
		  	gfPoly[] res1 = new gfPoly[res.length];
			for(int resLoc = 0; resLoc < res1.length; resLoc++)
			{
				res1[resLoc] = res[resLoc][0];
			}
			syndrome = new Poly(res1);
			   return syndrome;			/*
			Poly factor;
			gfPoly prod;
			Poly syndrome = new Poly(Poly.zeroPoly, 0);
			for(int i = 0; i < encoded.length; i++)
			{
				if(encoded[i]==1)
				{
					gfPoly toInv = Poly.eval(irredPoly, support[i], gF.gf_irredPoly);
					System.out.println("getting modular inverse of " + Poly.eval(irredPoly, support[i], gF.gf_irredPoly));
					prod = gfPoly.getModularInverse(Poly.eval(irredPoly, support[i], gF.gf_irredPoly), gF.gf_irredPoly);
					System.out.println("really inv? " + gfPoly.gf_multiply(toInv, prod, gF.gf_irredPoly));
					factor = Poly.getModularInverse(Poly.add(Poly.xPoly, new Poly(support[i])), irredPoly, gF.gf_irredPoly).get(1);
					factor = Poly.add(irredPoly, Poly.add(Poly.xPoly, new Poly(prod)));
					factor = Poly.multiply(factor, prod, gF.gf_irredPoly);
					factor = Poly.getRemainder(factor, irredPoly, gF.gf_irredPoly);
					syndrome = Poly.add(syndrome, factor);
				}
			}
			System.out.println("in here, got " + syndrome.toString());
			*/
			
		   return syndrome;
		}
		
		public static gfPoly[][] mult(int[] vec, Poly[][] polyPar, gfPoly fieldPol)
		{
			gfPoly[][] vecMat = new gfPoly[vec.length][1];
			for(int i = 0; i < vec.length; i++)
			{
			
				vecMat[i][0] = new gfPoly(vec[i]);
				
			}
			
			if(polyPar[0].length != vecMat.length)
			{
				return null; 
			}
			gfPoly[][] result = new gfPoly[polyPar.length][vecMat[0].length];
					
			
			
			int matrix1Row, matrix2Column, iterIndex;
			
			for(matrix1Row = 0; matrix1Row < polyPar.length; matrix1Row++)
			{
				for(matrix2Column = 0; matrix2Column < vecMat[0].length;
						matrix2Column++)
				{
					result[matrix1Row][matrix2Column] = new gfPoly(0);
					for(iterIndex = 0; iterIndex < vecMat.length;
							iterIndex++)
						{
							result[matrix1Row][matrix2Column] = gF.gf_add(result[matrix1Row][matrix2Column], Poly.multiply(polyParity[matrix1Row][iterIndex],
									new Poly(new gfPoly[] { vecMat[iterIndex][matrix2Column]}), gF.gf_irredPoly).coefficients[0]);
							

						}
					
				}
					
			}
			
			return result;
					
		}
		
}
