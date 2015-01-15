package utils;

public class SequenceUtils {
	public static String getRepresentative(String sequence){
		StringBuilder str_double = new StringBuilder();
		String temp;

		str_double.append(sequence);
		str_double.append(sequence);


		String representative = sequence;

		for (int i = 0; i < sequence.length(); i++) {
			temp = str_double.substring(i, i + sequence.length());
			if (temp.compareTo(representative) < 0) {
				representative = temp;
			}
		}

		str_double = str_double.reverse();

		temp = new String(str_double);

		temp = temp.replace('C', 'Z');
		temp = temp.replace('G', 'C');
		temp = temp.replace('Z', 'G');

		temp = temp.replace('T', 'Z');
		temp = temp.replace('A', 'T');
		temp = temp.replace('Z', 'A');

		str_double = new StringBuilder(temp);
		for (int i = 0; i < sequence.length(); i++) {
			temp = str_double.substring(i, i + sequence.length());
			if (temp.compareTo(representative) < 0) {
				representative = temp;
			}
		}


		return representative;
	}
	public static String getReverse(String sequence){
		StringBuilder str = new StringBuilder();
		str.append(sequence);
		return str.reverse().toString();
	}
	public static String getReverseComplementary(String sequence){
		String reversed = getReverse(sequence);
		StringBuilder builder = new StringBuilder();
		builder.append(reversed);
		
		String temp = new String(builder);
		temp = temp.replace('C', 'Z');
		temp = temp.replace('G', 'C');
		temp = temp.replace('Z', 'G');

		temp = temp.replace('T', 'Z');
		temp = temp.replace('A', 'T');
		temp = temp.replace('Z', 'A');
		
		return temp;
	}
}
