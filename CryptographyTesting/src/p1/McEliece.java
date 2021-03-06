package p1;

import java.lang.Math;
import java.math.BigInteger;
import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;

public class McEliece
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
   
   public static int[][] genMatrix;
   
   public static Poly[][] polyParity;
   
   public static int[] errorVec;
   
   public static int[][] pubKey;
   
   public static int[][] unreducedParity;
   
	public static int[][] createPublicKey(int numErrors, int extDegree, int supportSize)
	{
	   
	  
	   //irredPoly = Polynomial.getIrredPoly(numErrors); 
	   int[][] generatorMatrix = createGeneratorMatrix(numErrors, extDegree, supportSize);
	   genMatrix = generatorMatrix;
	   int[][] permMatrix = createPermMatrix(generatorMatrix[0].length);	 
	   
	   permutation = permMatrix;
	   
	   

		int rows = generatorMatrix.length;
		

		
		int[][] mA = new int[rows][rows];
		
		int[][] mT = new int[rows][rows];
		
	
		int[][] randomMatrix = getRandMatrix(mA, mT, rows);		
		
		scramble = randomMatrix;

		int[][] publicKey =  multiply(generatorMatrix, permMatrix);
				
		publicKey = multiply(randomMatrix, publicKey);
		
		return publicKey;
		

	}
	
	public static Boolean vecEquals(int[] one, int[] other)
	{
	   if(one.length != other.length)
	   {
	      return false;
	   }
	   for(int i = 0; i < other.length; i++)
	   {
	      if(one[i] != other[i])
	      {
	         return false;
	      }
	   }
	   return true;
	}
	
	public static int[] encrypt(int[] message, int numErrors, int extDegree, int suppSize )
	{
      
	   int[][] publicKey = createPublicKey(extDegree, numErrors, suppSize);
	   pubKey = publicKey;
	   

	   int[][] modMessage = new int[1][message.length];
	   modMessage[0] = message;

	   int[][] tempEnc = multiply(modMessage, publicKey);
	   
	   int[] error = getRandomErrorVector(numErrors, tempEnc[0].length);

	   
	   errorVec = error;
	
	   
	   
	   int[] encrypted = tempEnc[0];
	   for(int i = 0; i < encrypted.length; i++)
	   {
	      encrypted[i] ^= error[i];
	   }
	   
	   System.out.println("Encrypted:");
	   print(encrypted);
	   
	   return encrypted; //TODO:temp
	   
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
	
	public static int[][] createGeneratorMatrix(int extDegree, int numErrors, int supportSize)
	{
	   int[][] parityCheckMatrix;
	   do
	   {
	      parityCheckMatrix = getParityCheck(extDegree, numErrors, supportSize);
	   }while(!isSystematic(parityCheckMatrix, parityCheckMatrix[0].length-parityCheckMatrix.length));
	   parityCheck = parityCheckMatrix;
	   
	   
	   int[][] generatorMatrix = new int[parityCheckMatrix[0].length - parityCheckMatrix.length][support.length];
	   
	   
      for(int idRow = 0; idRow < generatorMatrix.length; idRow++)
      { 
         for(int idCol = 0; idCol < generatorMatrix.length; idCol++)
         {
            if(idRow == idCol)
            {
               generatorMatrix[idRow][idCol] = 1;
            }
         }
      }
	   for(int copyingCol = 0; copyingCol < parityCheckMatrix[0].length - parityCheckMatrix.length; copyingCol++)
	   {
	      for(int copyingRow = 0; copyingRow < parityCheckMatrix.length; copyingRow++)
	      { 
 	         generatorMatrix[copyingCol][copyingRow + generatorMatrix.length] 
	               = parityCheckMatrix[copyingRow][copyingCol];
	      }
	   }
	 
	         
      
      return generatorMatrix;
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
	
	
	public static int[] getRandomErrorVector(int weight, int length)
	{
		//error vector
		int[] vector = new int[length];
		
		//locations where 1s are placed
		int[] locations = new int[weight];

		for(int i = 0; i < weight; i++)
		{
			locations[i] = length;
		}
		
		int addTo, numAdded;
	
		for(numAdded = 0; numAdded < weight; numAdded++)
		{
			do
			{
		
				addTo = (int)(Math.random() * length);
		
			}while(isInArray(locations, addTo));
			
			locations[numAdded] = addTo;
			vector[addTo] = 1;
		}
		
		
		return vector;
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
	
   public static void print(gfPoly[][] array)
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
	public static gfPoly[][] multiply(gfPoly[][] matrix1, gfPoly[][] matrix2, gfPoly fieldPoly)
	{
		if(matrix1[0].length != matrix2.length)
		{
			return null; 
		}
		gfPoly[][] result = new gfPoly[matrix1.length][matrix2[0].length];
				
		
		int matrix1Row, matrix2Column, iterIndex;
		
		for(matrix1Row = 0; matrix1Row < matrix1.length; matrix1Row++)
		{
			for(matrix2Column = 0; matrix2Column < matrix2[0].length;
					matrix2Column++)
			{
				result[matrix1Row][matrix2Column] = new gfPoly(0);
				for(iterIndex = 0; iterIndex < matrix2.length;
						iterIndex++)
					{
					   result[matrix1Row][matrix2Column] = gfPoly.gf_add(result[matrix1Row][matrix2Column],
							gfPoly.gf_multiply(matrix1[matrix1Row][iterIndex], matrix2[iterIndex][matrix2Column], 
									fieldPoly)); 
		

					}

				
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
	
	public static int[][] getParityCheck(int extenDegree, int numErrors, int suppSize)
	{
       gF = new GaloisField(extenDegree);
       extDegree = extenDegree;
       McEliece.numErrors = numErrors;
      irredPoly = Poly.getIrredPoly(numErrors, gF.gf_irredPoly);
      support = gF.getSupport(irredPoly, suppSize);
      
      System.out.println("support:");
      print(support);

      
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
	/*
	public static List<Polynomial> berlTrace(Polynomial toFactor, Polynomial galoisMod)
	{
	   int r = 1;
	   Polynomial h = new Polynomial();
	   Polynomial trace = getTracePolynomial();
	   
	   
	   List<Polynomial> factors = new ArrayList<Polynomial>();
	   factors.add(toFactor);
	   
	   for(int j = 0; j < extDegree; j++)
	   {
	      if(r != toFactor.degree)
	      {
	         h = Polynomial.eval(trace, inner, mod)
	      }
	   }
	}
	*/
	
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
	
	public static Poly getSuppInvTest(int suppLoc)
	{
	   Poly inv = Poly.zeroPoly;
	   Poly fact = new Poly();
	   
	   
	   for(int parityLoc = 0; parityLoc < polyParity.length; parityLoc++)
	   {
         fact = Poly.toPower(Poly.xPoly, gF.gf_irredPoly, polyParity.length - parityLoc - 1);

         fact = Poly.multiply(fact, polyParity[parityLoc][suppLoc], gF.gf_irredPoly);
         inv = Poly.add(inv, fact);
         
	   }
	   
	   
	   Poly poly = Poly.xPoly;
	   Poly supp = new Poly(support[suppLoc]);
	   poly = poly.add(poly, supp);
	   
	   return inv;
	}
	
	public static Poly altSyndrome(int[] encoded)
	{
	   Poly factor;
	   Poly syndrome = new Poly(Poly.zeroPoly, 0);
	   for(int i = 0; i < encoded.length; i++)
	   {
	      if(encoded[i] == 1)
	      {
	    	  Poly toInv = Poly.add(Poly.xPoly, new Poly(support[i]));
System.out.println("at i= " + i + " getting inverse of " + toInv.toString());
Poly inv = Poly.getModularInverse(toInv, irredPoly, gF.gf_irredPoly).get(1);
System.out.println("inverse is: " + inv.toString());
System.out.println("got inv? " + Poly.multiply(toInv, inv, gF.gf_irredPoly, irredPoly).toString());
	         factor = Poly.getModularInverse(Poly.add(Poly.xPoly, new Poly(support[i])), irredPoly, gF.gf_irredPoly).get(1);
	         syndrome = Poly.add(syndrome,  factor);
	      }
	      
	      
	   }
	   
	   return syndrome;
	}
	
	
	
}

 
