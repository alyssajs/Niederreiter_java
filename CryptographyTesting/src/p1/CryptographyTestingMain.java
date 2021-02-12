package p1;

public class CryptographyTestingMain {
	
	public static void main(String[] args)
	{
		int[] msg = {0, 1, 0};
				
		Message test = new Message(msg);
						
		int[] encoded = test.mapToCodeWord(5, 3);
		
		/*int index;
		for(index = 0; index < encoded.length; index++)
		{
			System.out.print(encoded[index] + " " );
		}
		*/
		
		/*int[][] randomMatrix = createRandomMatrixS(5);
		
		for(int i = 0; i < 5; i++)
		{
			for(int j = 0; j < 5; j++)
			{
				System.out.print(randomMatrix[i][j] + " " );
			}
			
			System.out.println();
		}
		
		int[] errorVector = getRandomErrorVector(3, 5);
		System.out.println("Error vector:");
		for(int i = 0; i < errorVector.length; i++)
		{
			System.out.print(errorVector[i] + " ");
		}*/
		
		int[][] permMatrix = createPermMatrix(5);
		
		System.out.println("Permutation matrix:");
		for(int i = 0; i < 5; i++)
		{
			for(int j = 0; j < 5; j++)
			{
				System.out.print(permMatrix[i][j] + " ");
			}
			System.out.println();
		}
		
		
	}
	
	public static int[][] createRandomMatrixS(int dim)
	{
		int[][] randomMatrix = new int[dim][dim];
		
		int row, column, entry;
		
		for(row = 0; row < dim; row++)
		{
			for(column = 0; column < dim; column++)
			{
				entry = (int)(Math.random() * 2);
				randomMatrix[row][column] = entry;
			}
		}
		
		return randomMatrix;
	}
	
	public static int[] getRandomErrorVector(int weight, int length)
	{
		int[] vector = new int[length];
		
		int[] locations = {length, length, length};
		
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
	
	public static int[][] createPermMatrix(int length)
	{
		int[][] permMatrix = new int[length][length];
		
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

}
