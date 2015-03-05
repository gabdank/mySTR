package domain;

import org.biojava3.alignment.template.SequencePair;
import org.biojava3.core.sequence.DNASequence;
import org.biojava3.core.sequence.compound.NucleotideCompound;

public class AlignmentsData {

	private SequencePair<DNASequence, NucleotideCompound> leftSideAlignment;
	private SequencePair<DNASequence, NucleotideCompound> rightSideAlignment;
	private SequencePair<DNASequence, NucleotideCompound> pairEndAlignment;
	
	private String leftSeq;
	private String rightSeq;
	private String repetitiveSequence;
	private String readID;
	private double score;

	public AlignmentsData(double scr, String id, String l, String r, String repeat,SequencePair<DNASequence, NucleotideCompound> left,
			SequencePair<DNASequence, NucleotideCompound> right, SequencePair<DNASequence, NucleotideCompound> pair ){  
		leftSeq = l;
		rightSeq = r;
		leftSideAlignment = left;
		rightSideAlignment = right;
		repetitiveSequence =repeat;
		pairEndAlignment = pair;
		readID = id;
		score = scr;
	}

	public String getLeftSeq(){
		return leftSeq;
	}

	public String getRigthSeq(){
		return rightSeq;
	}

	public SequencePair<DNASequence, NucleotideCompound> getLeftSideAlignment() {
		return leftSideAlignment;
	}

	public SequencePair<DNASequence, NucleotideCompound> getRightSideAlignment() {
		return rightSideAlignment;
	}
	
	public SequencePair<DNASequence, NucleotideCompound> getPairEndAlignment() {
		return pairEndAlignment;
	}
	
	public String getRepetitiveSequence() {
		return repetitiveSequence;
	}

	public String getReadID() {
		return readID;
	}

	public double getScore() {
		return score;
	}

	public void setLeftSide(SequencePair<DNASequence, NucleotideCompound> leftSide) {
		this.leftSideAlignment = leftSide;
	}

	public void setRightSide(SequencePair<DNASequence, NucleotideCompound> rightSide) {
		this.rightSideAlignment = rightSide;
	}

	public void setRepetitiveSequence(String repetitiveSequence) {
		this.repetitiveSequence = repetitiveSequence;
	}

	public void setReadID(String readID) {
		this.readID = readID;
	}

	public void setScore(double score) {
		this.score = score;
	}

	public String toString(){
		return "Read ID:"+ readID+"\t ALIGNMENT SCORE:"+score+"\nLEFT SIDE:\n"+leftSideAlignment+"\n"+leftSeq
				+"\nRIGHT SIDE\n"+rightSideAlignment+"\n"+rightSeq +"\nREPETITIVE SEQUECE\n"+repetitiveSequence;
	}

	public String getStringRep(int repetitiveUnitLength, int leftFlankMax, int rightFlankMax) {

		int numOfRepetitions = repetitiveSequence.length()/repetitiveUnitLength;

		String leftToPrint="";
		if (leftSeq.length()<leftFlankMax){
			int delta = leftFlankMax-leftSeq.length();
			for (int s = 0;s<delta;s++){
				leftToPrint += "_";
			}
			leftToPrint+=leftSeq;
		}
		else{
			leftToPrint=leftSeq.substring(leftSeq.length()-leftFlankMax);
		}

		String rightToPrint="";
		if (rightSeq.length()<rightFlankMax){
			int delta = rightFlankMax-rightSeq.length();
			rightToPrint+=rightSeq;
			for (int s = 0;s<delta;s++){
				rightToPrint += "_";
			}

		}
		else{
			rightToPrint=rightSeq.substring(0, rightFlankMax);
		}

		String repetitiveUnit = repetitiveSequence.substring(0, repetitiveUnitLength);


		String answer = leftToPrint+"\t{"+repetitiveUnit+"}"+numOfRepetitions+"\t"+rightToPrint;
		return answer;
	}
	
	

	public String reportTheAlignmentData(){
		String toReturn = "";
		toReturn = "leftSequence:"+leftSeq+"$"+"rightSequence:"+rightSeq+"$"+"repetitiveSequence:"+repetitiveSequence+"$"+"readID:"+readID+"$"+
		"alignmentScore:"+score+"\n";
		toReturn+="leftAlignment:"+leftSideAlignment.getQuery()+"%"+leftSideAlignment.getTarget()+"$rightAlignment:"+rightSideAlignment.getQuery()+
				"%"+rightSideAlignment.getTarget()+"$pairedEndAlignment:"+pairEndAlignment.getQuery()+"%"+pairEndAlignment.getTarget();
		
		
		return toReturn;
	}
	

}
