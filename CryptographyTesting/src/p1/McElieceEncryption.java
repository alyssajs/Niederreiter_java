package p1;

import java.lang.Math;
import java.util.Scanner;

public class McElieceEncryption {
/*
	public static void main(String[] args)
	{
		int[] msg = {0,1,1,0};
		
		genMatrix = createGeneratorMatrix)(4, 8);
		
	}
	
	public int[][] createGeneratorMatrix(int rows, int columns)
	{
		int rowIndex; columnIndex;
		
		Scanner scan = new Scanner();
		
		int[][] generatorMatrix = new int[rows][columns];
		
		for(rowIndex = 0; rowIndex < rows; rowIndex++)
		{
			for(columnIndex = 0; columnIndex < columns; columnIndex++)
			{
				generatorMatrix[rowIndex][columnIndex] = scan.nextInt();
			}
		}
		
		return generatorMatrix;
	}

	public int[][] createRandomMatrixS(int dim)
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
	
	
	public int[] getRandomErrorVector(int weight, int length)
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
		
				addTo = (int)(Math.random() * length)
		
			}while(isInArray(locations, addTo))
			
			vector[addTo] = 1;
		}
		
		
		return vector;
	}
	
	public boolean isInArray(int[] array, int testVal)
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
	
	public int[][] getPublicKey()
	{
		//return S * G * P
	}
	
	public int[][] createIDMatrix(int length)
	{
		int index;
		int[][] idMatrix = new int[][];
		
		for(index = 0; index < length; index++)
		{
		
			idMatrix[index][index] = 1;
		
		}
		
		return idMatrix;
	}
	
	public int[] encryptMcEliece()
	{
		
	}
	*/
}

