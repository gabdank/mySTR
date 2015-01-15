package domain.repeatDetectionMorishita;

public class ThinnerStackRecord {
	private final int startIndex;
	private final int endIndex;
	private final int treeDepth;
	private final int treeDepthZero;
	
	public ThinnerStackRecord(int currentStart, int currentEnd, int depth, int depthZero){
		startIndex = currentStart;
		endIndex = currentEnd;
		treeDepth = depth;
		treeDepthZero = depthZero;
		
	}

	

	public int getStartIndex() {
		return startIndex;
	}
	
	public int getEndIndex() {
		return endIndex;
	}
	
	public int getTreeDepth() {
		return treeDepth;
	}

	public int getTreeDepthZero() {
		return treeDepthZero;
	}
	
	public String toString(){
		return startIndex + "\t" +endIndex+"\t"+ treeDepth + "\t" + treeDepthZero;  
	}
}
