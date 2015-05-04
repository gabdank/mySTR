package utils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.StringTokenizer;

import domain.repeatDetectionMorishita.FindRepetition;
import domain.repeatDetectionMorishita.RepetitionList;

public class bedFileFilter {

	public static void main(String[] args) throws IOException {

		BufferedReader  br = new BufferedReader(new FileReader("/home/gabdank/Documents/ElegansSTRValidation/WS243.two.bed"));
		//BufferedReader  br = new BufferedReader(new FileReader("trfIndexNEW1.bed"));
		String line;
		int mone =0;

		System.out.println("STARTED READING BED FILE");


		FileWriter fw1 = new FileWriter(new File("/media/gabdank/Disk3/mySTR/WS243.two.filtered.bed"));
		BufferedWriter bw1 = new BufferedWriter(fw1);


		while ((line=br.readLine())!=null){
			StringTokenizer t = new StringTokenizer(line);
			int numberOfTokens = t.countTokens();
			//System.out.println("Number of tokens: "+ numberOfTokens);
			if (numberOfTokens == 7){
				String chromosome = t.nextToken();
				int start = Integer.parseInt(t.nextToken());
				int end = Integer.parseInt(t.nextToken());

				String unit=t.nextToken();
				String sequence=t.nextToken();
				int repLength = Integer.parseInt(t.nextToken());
				int startInside = Integer.parseInt(t.nextToken());
				
				String representativeUnit = SequenceUtils.getRepresentative(unit);

				
				int repeatLength = end-start+1;
				
				System.out.println(line);
				
				String STR_Region = sequence.substring(startInside-1, startInside+repeatLength-1);
				String leftFlank = sequence.substring(0, startInside-1);
				String rightFlank = sequence.substring(startInside+repeatLength-1);
				System.out.println("-----------------");
				System.out.println(STR_Region);				
				System.out.println(leftFlank);				
				System.out.println(rightFlank);
				System.out.println("-----------------");
				
				if (leftFlank.length()>20){
					leftFlank = leftFlank.substring(leftFlank.length()-20);
					System.out.println(leftFlank);
				}
				if (rightFlank.length()>20){
					rightFlank = rightFlank.substring(0,20);
					System.out.println(rightFlank);
				}
				System.out.println("================================");
				String combinedRepeat = leftFlank+STR_Region+rightFlank;
				RepetitionList repetition1 = new RepetitionList();
				FindRepetition find_repetition= new FindRepetition();
				repetition1.emptyRepetitionList();
				find_repetition.findMaximalRepetition(combinedRepeat,repetition1);
				find_repetition.extendMaximalRepetition("sample_id", combinedRepeat, repetition1);
				
				//System.out.println("unit size : "+repetition1.get_max_repetition_unit_size());
				//System.out.println("length of the repetition : "+repetition1.get_max_repetition_length());
				// Ensure we are dealing with short tandem repeat (2..8bp long) and that the repeat length is at least 16bp long 
				if (repetition1.get_max_repetition_unit_size()>1 && repetition1.get_max_repetition_unit_size()<9 && repetition1.get_max_repetition_length()>16){
					String repeat1 = combinedRepeat.substring(repetition1.get_max_repetition_left_index(),repetition1.get_max_repetition_left_index()+repetition1.get_max_repetition_unit_size());
					System.out.println("Repeat = "+repeat1);
				}
				System.out.println("================================");
					
					
				
				/** 
				 * What do we have to check?:
				 * 1) connect flanks and run through the repeat detector.
				 * 
				 */
				
			if (mone==3)
				break;
			mone ++;
				
			//if (mone%100000==0) 
			//	System.out.println("PROCESSED "+mone +" lines in BED");
			//mone++;
			}
		}

	}

}
