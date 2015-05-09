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

		BufferedReader  br = new BufferedReader(new FileReader("/home/gabdank/Documents/HumanSTR/hg19.two.bed"));
		//BufferedReader  br = new BufferedReader(new FileReader("trfIndexNEW1.bed"));
		String line;
		int mone =0;

		System.out.println("STARTED READING BED FILE");


		FileWriter fw1 = new FileWriter(new File("/home/gabdank/Documents/HumanSTR/hg19.two.filtered.bed"));
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

				
				//int repeatLength = end-start+1;
				//System.out.println(repLength + "?=?"+repeatLength);
				
				//System.out.println(line);
				/*if (startInside+repLength==769){
					System.out.println(line);
					System.out.println("sequence length: "+sequence.length());
					System.out.println(startInside+repLength);
					System.out.println(sequence.substring(startInside-1, startInside+repLength-1));
					System.out.println((startInside-1)+"\t"+(startInside+repLength-1)+"\t"+sequence.length());
				}*/
				String STR_Region = sequence.substring(startInside-1, startInside+repLength-1);
				String leftFlank = sequence.substring(0, startInside-1);
				String rightFlank = sequence.substring(startInside+repLength-2);
				
				
				if (leftFlank.length()>30){
					leftFlank = leftFlank.substring(leftFlank.length()-30);
					//System.out.println(leftFlank);
				}
				if (rightFlank.length()>30){
					rightFlank = rightFlank.substring(0,30);
					//System.out.println(rightFlank);
				}
				// Creating the combined sequence - that contains flanks and the STR region between them
				String combinedRepeat = leftFlank+STR_Region+rightFlank;
				int startPosition = leftFlank.length();
				int endPosition = startPosition+STR_Region.length();
				
				//System.out.println("REFERENCE:");
				//System.out.println(leftFlank+">>"+STR_Region+"<<"+rightFlank);
				
				//System.out.println("Updated Start = "+startPosition+ " Updated  End = "+ endPosition);
				
				
				
				
				RepetitionList repetition1 = new RepetitionList();
				FindRepetition find_repetition= new FindRepetition();
				repetition1.emptyRepetitionList();
				find_repetition.findMaximalRepetition(combinedRepeat,repetition1);
				find_repetition.extendMaximalRepetition("sample_id", combinedRepeat, repetition1);
				if (repetition1.get_max_repetition_unit_size()>1 && repetition1.get_max_repetition_unit_size()<9 && repetition1.get_max_repetition_length()>16){
					String repeat1 = combinedRepeat.substring(repetition1.get_max_repetition_left_index(),repetition1.get_max_repetition_left_index()+repetition1.get_max_repetition_unit_size());
					String repeatSection = combinedRepeat.substring(repetition1.get_max_repetition_left_index(),repetition1.get_max_repetition_right_index() );
					
					int startPositionNew = repetition1.get_max_repetition_left_index();
					int endPositionNew = repetition1.get_max_repetition_right_index();
					
					//System.out.println("New Updated Start = "+startPositionNew+ " New Updated  End = "+ endPositionNew);
					
					int startDelta = startPositionNew-startPosition; // check if it shouldn't be the other way around
					int endDelta = endPositionNew-endPosition;
					
					//System.out.println("Start delta = "+startDelta+" End delta = "+endDelta);
					
					//System.out.println(combinedRepeat.substring(0,repetition1.get_max_repetition_left_index())+"**"+repeatSection+"**"+combinedRepeat.substring(repetition1.get_max_repetition_right_index()));
					
					//System.out.println(chromosome+"\t"+(start+startDelta) + "\t" +(end+endDelta));
										
					bw1.write(chromosome+"\t"+(start+startDelta)+"\t"+(end+endDelta)+"\t"+repeat1+"\t"+sequence+"\t"+repeatSection.length()+"\t"+(startInside+startDelta)+"\n");
					
				}
				/*else{
					System.out.println("WE HAVE A PROBLEM");
					System.out.println(line);
					System.out.println(repetition1.get_max_repetition_length());
					System.out.println(repetition1.get_max_repetition_unit_size());
					String repeat1 = combinedRepeat.substring(repetition1.get_max_repetition_left_index(),repetition1.get_max_repetition_left_index()+repetition1.get_max_repetition_unit_size());
					System.out.println(repeat1);
					throw new IOException("ZOPA");
				}*/
			
				//System.out.println("================================");
					
					
				
				/** 
				 * What do we have to check?:
				 * 1) connect flanks and run through the repeat detector.
				 * 
				 */
				
			//if (mone==3)
			//	break;
			mone ++;
				
			if (mone%100000==0) 
				System.out.println("PROCESSED "+mone +" lines in BED");
			mone++;
			}
		}
		bw1.close();
		fw1.close();

	}

}
