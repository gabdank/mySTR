package domain.repeatDetectionMorishita;



import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Iterator;

public class RepetitionList {
	 private ArrayList<int[]> List = new ArrayList<int[]>();
		
		private int maxRepetitionLength = 0;
		private int maxRepetitionUnitLength = 0;
		private int maxRepetitionUnitLeft =0;
		private int maxRepetitionUnitRight = 0;; 

	public void emptyRepetitionList(){
		 List = new ArrayList<int[]>();
		 maxRepetitionLength = 0;
		 maxRepetitionUnitLength = 0;
		 maxRepetitionUnitLeft =0;
	     maxRepetitionUnitRight = 0; 

	}

   
	public void addRepetition(int j, int left, int right) {
		for (int i = 0; i < List.size(); i++) {
		    int[] ob = List.get(i);
		    if (ob[1] == left && ob[2] == right) {
			return;
		    }
		}
		int[] ob = {j,left,right};
		if (maxRepetitionLength < right - left + 1){
			maxRepetitionLength = right -left +1;
			maxRepetitionUnitLength = j;
			maxRepetitionUnitLeft = left;
			maxRepetitionUnitRight = right;
		}
		List.add(ob);
	}

	public void print_repetition() {
		for (int i = 0; i < List.size(); i++) {
			System.out.println("Period : " + List.get(i)[0] + " Left : "
					+ List.get(i)[1] + " Right : " + List.get(i)[2]);
		}
	}
	
	public void print_repetition(PrintWriter file) {
		for (int i = 0; i < List.size(); i++) {
			file.println("Period : " + List.get(i)[0] + " Left : "
					+ List.get(i)[1] + " Right : " + List.get(i)[2]);
		}
	}

	public void print_repetition(String text) {
		for (int i = 0; i < List.size(); i++) {
		    int[] ob = List.get(i);
		    System.out.println("Period : "
				       + ob[0]
				       + " Left : "
				       + ob[1]
				       + " Right : "
				       + ob[2]
				       + " repetition : "
				       + text.substring(ob[1], ob[1] + ob[0]));
		}
	}
	
	public void print_repetition(PrintWriter file, String text) {
		for (int i = 0; i < List.size(); i++) {
		    int[] ob = List.get(i);
		    file.println("Period : "
				       + ob[0]
				       + " Left : "
				       + ob[1]
				       + " Right : "
				       + ob[2]
				       + " repetition : "
				       + text.substring(ob[1], ob[1] + ob[0]));
		}
	}

	public String get_repetition_string(String text, int index) {
	    int[] ob = List.get(index);
	    return text.substring(ob[1], ob[1] + ob[0]);
	}

	public int get_MaximalRepetition_length(int index) {
	    int[] ob = List.get(index);
		return ob[2] - ob[1] + 1;
	}

	public int get_Repetition_ListSize() {
		return List.size();
	}

        public int  get_repetition_left_index(int index) {
	    return List.get(index)[1];
	}

        public int  get_repetition_right_index(int index) {
	    return List.get(index)[2];
	}

        public int  get_repetition_unit_size(int index) {
	    return List.get(index)[0];
	}
        
    public int  get_max_repetition_length() {
        	return maxRepetitionLength;
	}
    public int  get_max_repetition_left_index() {
        	return maxRepetitionUnitLeft;
	}

    public int  get_max_repetition_right_index() {
        	return maxRepetitionUnitRight;
	}

    public int  get_max_repetition_unit_size() {
        	return maxRepetitionUnitLength;
	}
        
        
    public void insert(int j, int left, int right) {
    	Iterator<int[]> itr = List.iterator();
    	while(itr.hasNext()){
    		int[] ob = (int[]) itr.next();
    		if(ob[0] == j && ob[1] >= left && ob[2] <= right){
    			itr.remove();
    		}
    	}
    	int[] o = {j,left,right};
    	List.add(o);
		maxRepetitionLength = right -left +1;
		maxRepetitionUnitLength = j;
		maxRepetitionUnitLeft = left;
		maxRepetitionUnitRight = right;
    }

}
