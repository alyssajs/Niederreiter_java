package p1;

public class Message {

	public int[] message;
	
	public Message(int[] msg)
	{
		message = msg;
	}

	public int[] encryptNiederetter(int length, int w)
	{
		int[] encrypted = new int[length];
		
		return encrypted;
	}
	
	/**
	 * Maps a bit string to a constant weight codeword
	 * @param length length of output
	 * @param weight weight of code word
	 * @return encoded message
	 */
	public int[] mapToCodeWord(int length, int weight)
	{
		int[] codeWord = new int[length];
		
		double c, cPrime, i;
		
		int j;
		
		i = calcLexIndex();
		
		c = factorial(length) / 
				(factorial(weight) * factorial(length - weight));
		
		cPrime  = 0;
		
		j = length;
		
		while(j > 0)
		{
			cPrime = c * (j - weight) / j;
			
			if(i <= cPrime)
			{
				codeWord[ j - 1 ] = 0;
				
				c = cPrime;
			}
			
			else
			{
				codeWord[ j - 1 ] = 1;
				
				i = i - cPrime;
				
				c = c * weight / length;
		
			}
			
			j--;
		}
		
		return codeWord;
	}
	
	public int factorial(int input)
	{
		if(input > 1)
		{
			return input * factorial(input - 1);
		}
		return 1;
	}
	
	public int calcLexIndex()
	{
		int binaryIndex, result;
		
		result = 0;
		
		for(binaryIndex = message.length - 1; binaryIndex >= 0; binaryIndex--)
		{
			result += message[binaryIndex] 
						* Math.pow(2, message.length - binaryIndex - 1);
		}
		
		return result + 1;
	}
	
	
}
