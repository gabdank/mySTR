package domain.repeatDetectionMorishita;



import java.util.Random;


public class FindRepetition {
	private int[] calc_LP(String text, boolean isReverse) {
		int length = text.length();
		char[] v;
		int[] lp;
		if (length % 2 == 0) {
			lp = new int[length / 2];
			v = new char[length / 2];
			for (int i = 0; i < length / 2; i++) {
				v[i] = text.charAt(i + length / 2);
			}
		} else {
			if (!isReverse) {
				lp = new int[(length + 1) / 2];
				v = new char[(length + 1) / 2];
				for (int i = 0; i < (length + 1) / 2; i++) {
					v[i] = text.charAt(i + length / 2);
				}
			} else {
				lp = new int[(length - 1) / 2];
				v = new char[(length - 1) / 2];
				for (int i = 0; i < (length - 1) / 2; i++) {
					v[i] = text.charAt(i + length / 2 + 1);
				}
			}
		}
		lp = calc_lppattern(String.valueOf(v));
		return lp;
	}

	private int[] calc_LS(String text, boolean isReverse) {
		int length = text.length();
		int[] temp;

		if (!isReverse) {
			char[] u = new char[length / 2];
			for (int i = 0; i < length / 2; i++) {
				u[i] = text.charAt(i);
			}

			// System.out.println("u : " + String.valueOf(u));

			char[] inverseU = new char[length / 2];
			for (int i = 0; i < length / 2; i++) {
				inverseU[i] = u[u.length - i - 1];
			}

			char[] inverseUV = new char[text.length()];
			for (int i = 0; i < text.length(); i++) {
				inverseUV[i] = text.charAt(text.length() - i - 1);
			}

			int[] lppattern = calc_lppattern(String.valueOf(inverseU));
			int[] lptext = calc_lptext(String.valueOf(inverseUV),
					String.valueOf(inverseU), lppattern);

			// System.out.println("text : " + String.valueOf(inverseUV));
			// System.out.println("patr : " + String.valueOf(inverseU));

			int[] ls = new int[text.length() - text.length() / 2];
			for (int i = 0; i < text.length() - text.length() / 2; i++) {
				ls[i] = lptext[i];
			}

			temp = new int[ls.length];
			for (int i = 0; i < ls.length; i++) {
				temp[i] = ls[ls.length - i - 1];
			}
		} else {

			char[] u = new char[length / 2 + length % 2];
			for (int i = 0; i < length / 2 + length % 2; i++) {
				u[i] = text.charAt(i);
			}

			// System.out.println("u : " + String.valueOf(u));

			char[] inverseU = new char[length / 2 + length % 2];
			for (int i = 0; i < length / 2 + length % 2; i++) {
				inverseU[i] = u[u.length - i - 1];
			}

			char[] inverseUV = new char[text.length()];
			for (int i = 0; i < text.length(); i++) {
				inverseUV[i] = text.charAt(text.length() - i - 1);
			}

			int[] lppattern = calc_lppattern(String.valueOf(inverseU));
			int[] lptext = calc_lptext(String.valueOf(inverseUV),
					String.valueOf(inverseU), lppattern);

			// System.out.println("text : " + String.valueOf(inverseUV));
			// System.out.println("patr : " + String.valueOf(inverseU));

			int[] ls = new int[text.length() - text.length() / 2 - length % 2];
			for (int i = 0; i < ls.length; i++) {
				ls[i] = lptext[i];
			}

			temp = new int[ls.length];
			for (int i = 0; i < ls.length; i++) {
				temp[i] = ls[ls.length - i - 1];
			}

		}

		return temp;
	}

	private int[] calc_lptext(String text, String pattern, int[] lppattern) {

		int[] lptext = new int[text.length()];

		// First find lppattern[1]
		int j = 0;
		if (j < text.length() && j < pattern.length()) {
			while (text.charAt(j) == pattern.charAt(j)) {
				j++;
				if (j == text.length() || j == pattern.length()) {
					break;
				}
			}
		}
		lptext[0] = j;

		for (int i = 1; i < lptext.length; i++) {
			// calculate max length
			int k = 0;
			int max_length = 0;
			for (int l = 0; l < i; l++) {
				int length = l + lptext[l];
				if (max_length < length) {
					max_length = length;
					k = l;
				}
			}

			if (i < k + lptext[k]) {
				int x = max_length - i;

				if (x > lppattern[i - k]) {
					lptext[i] = lppattern[i - k];
				} else {
					j = 0;
					if ((x + i + j) < text.length()
							&& (x + j) < pattern.length()) {
						while (text.charAt(x + i + j) == pattern.charAt(x + j)) {
							j++;
							if (x + i + j == text.length()
									|| x + j == pattern.length()) {
								break;
							}
						}
					}
					lptext[i] = j + x;
				}
			} else {
				j = 0;
				while (pattern.charAt(j) == text.charAt(i + j)) {
					j++;
					if (i + j == text.length() || j == pattern.length()) {
						break;
					}
				}
				lptext[i] = j;

			}
		}

		return lptext;
	}

	private int[] calc_lppattern(String pattern) {

		int[] lppattern = new int[pattern.length()];

		lppattern[0] = pattern.length();
		if (pattern.length() == 1) {
			return lppattern;
		}
		// First find lppattern[1]
		int j = 0;
		if (j + 1 < pattern.length()) {
			while (pattern.charAt(j) == pattern.charAt(j + 1)) {
				j++;
				if (j + 1 == pattern.length()) {
					break;
				}
			}
		}
		lppattern[1] = j;

		for (int i = 2; i < lppattern.length; i++) {
			// calculate max length
			int k = 0;
			int max_length = 0;
			for (int l = 1; l < i; l++) {
				int length = l + lppattern[l];
				if (max_length < length) {
					max_length = length;
					k = l;
				}
			}

			if (i < k + lppattern[k]) {
				int x = max_length - i;

				if (x > lppattern[i - k]) {
					lppattern[i] = lppattern[i - k];
				} else {
					j = 0;
					if (x + i + j < pattern.length()) {
						while (pattern.charAt(x + i + j) == pattern.charAt(x
								+ j)) {
							j++;
							if (x + i + j == pattern.length()) {
								break;
							}
						}
					}

					lppattern[i] = j + x;
				}
			} else {
				j = 0;
				while (pattern.charAt(j) == pattern.charAt(i + j)) {
					j++;
					if (i + j == pattern.length()) {
						break;
					}
				}
				lppattern[i] = j;

			}
		}

		return lppattern;
	}

	public void findMaximalRepetition(String text, RepetitionList repetition_list) {
		findMaximalRepetition_rec(text, 0, 0, repetition_list, 0);
	}

	private void findMaximalRepetition_rec(String text, int left_index,	int tree_depth, RepetitionList repetition_list,
			int depth0_text_length) {
		//System.out.println("TEXT = "+text);
		//System.out.println("depth0 = "+depth0_text_length);
		if (text.length() == 1) {
			return;
		}
		char[] reverseText = new char[text.length()];
		for (int i = 0; i < text.length(); i++) {
			reverseText[i] = text.charAt(text.length() - i - 1);
		}

		// System.out.println("text : " + text);

		String u;
		if (text.length() > 1) {
			u = text.substring(0, text.length() / 2);
		} else {
			u = "";
		}
		String v = text.substring(text.length() / 2);

		// System.out.println("u : " + String.valueOf(u));
		// System.out.println("v : " + String.valueOf(v));

		int[] LP_temp = calc_LP(text, false);
		int[] LS_temp = calc_LS(text, false);

		int[] LP = new int[LP_temp.length + 1];
		int[] LS = new int[LS_temp.length + 1];

		for (int i = 0; i < LP_temp.length; i++) {
			LP[i] = LP_temp[i];
		}
		LP[LP.length - 1] = 0;

		for (int i = 0; i < LS_temp.length; i++) {
			LS[i] = LS_temp[i];
		}
		LS[LS.length - 1] = 0;

		for (int j = 1; j <= text.length() - text.length() / 2; j++) {
			if (LS[j - 1] + LP[j] >= j) {
				int left = text.length() / 2 - LS[j - 1];
				int right = text.length() / 2 + j + LP[j] - 1;

				if (tree_depth == 0
						|| !(left == 0 && left_index != 0)
						&& !(right == text.length() - 1 && (left_index + text
								.length()) != depth0_text_length)) {
					left += left_index;
					right += left_index;
					repetition_list.addRepetition(j, left, right);
				}
			}
		}
		int[] reverseLP = calc_LP(String.valueOf(reverseText), true);
		int[] reverseLS = calc_LS(String.valueOf(reverseText), true);


		for (int j = 1; j < text.length() - text.length() / 2 - text.length()
				% 2; j++) {


			if (reverseLS[j - 1] + reverseLP[j] >= j /* && reverseLS[j]<j */) {
				if(reverseLS[j-1] >= j){
					continue;
				}

				int right = reverseText.length
						- (reverseText.length / 2 - reverseLS[j - 1]) - 1
						- text.length() % 2;
				int left = reverseText.length
						- (reverseText.length / 2 + j + reverseLP[j] - 1)
						- 1 - text.length() % 2;

				if (tree_depth == 0
						|| !(left == 0 && left_index != 0)
						&& !(right == text.length() - 1 && (left_index + text
								.length()) != depth0_text_length)) {
					left += left_index;
					right += left_index;

					repetition_list.addRepetition(j, left, right);
				}

			}

		}
		if (tree_depth == 0) {
			findMaximalRepetition_rec(u, left_index, tree_depth + 1,
					repetition_list, text.length());
		} else {
			findMaximalRepetition_rec(u, left_index, tree_depth + 1,
					repetition_list, depth0_text_length);
		}

		if (tree_depth == 0) {
			findMaximalRepetition_rec(v, left_index + u.length(),
					tree_depth + 1, repetition_list, text.length());
		} else {
			findMaximalRepetition_rec(v, left_index + u.length(),
					tree_depth + 1, repetition_list, depth0_text_length);
		}
	}


/*	private void findMaximalRepetition_iter(String text, int left_index, int tree_depth, RepetitionList repetition_list,
			int depth0_text_length) {

		Stack localStack = new Stack();
		StackRecord record = new StackRecord(text, left_index, tree_depth, depth0_text_length);
		localStack.push(record);

		while (!localStack.isEmpty()){

			StackRecord currentRecord = localStack.pop();
			String currentText = currentRecord.getText();
			if (currentText.length()>1){
				int currentLeft = currentRecord.getLeftIndex();
				int currentTreeDepth = currentRecord.getTreeDepth();
				int currentZeroDepth = currentRecord.getTreeDepthZero();

				char[] reverseText = new char[currentText.length()];
				for (int i = 0; i < currentText.length(); i++) {
					reverseText[i] = currentText.charAt(currentText.length() - i - 1);
				}


				String u;
				if (currentText.length() > 1) {
					u = currentText.substring(0, currentText.length() / 2);
				} else {
					u = "";
				}
				String v = currentText.substring(currentText.length() / 2);



				int[] LP_temp = calc_LP(currentText, false);
				int[] LS_temp = calc_LS(currentText, false);

				int[] LP = new int[LP_temp.length + 1];
				int[] LS = new int[LS_temp.length + 1];

				for (int i = 0; i < LP_temp.length; i++) {
					LP[i] = LP_temp[i];
				}
				LP[LP.length - 1] = 0;

				for (int i = 0; i < LS_temp.length; i++) {
					LS[i] = LS_temp[i];
				}
				LS[LS.length - 1] = 0;

				for (int j = 1; j <= currentText.length() - currentText.length() / 2; j++) {
					if (LS[j - 1] + LP[j] >= j) {
						int left = currentText.length() / 2 - LS[j - 1];
						int right = currentText.length() / 2 + j + LP[j] - 1;

						if (currentTreeDepth == 0
								|| !(left == 0 && currentLeft != 0)
								&& !(right == text.length() - 1 && (currentLeft + text
										.length()) != currentZeroDepth)) {
							left += currentLeft;
							right += currentLeft;
							repetition_list.addRepetition(j, left, right);
						}
					}
				}
				int[] reverseLP = calc_LP(String.valueOf(reverseText), true);
				int[] reverseLS = calc_LS(String.valueOf(reverseText), true);


				for (int j = 1; j < currentText.length() - currentText.length() / 2 - currentText.length()
						% 2; j++) {


					if (reverseLS[j - 1] + reverseLP[j] >= j ) {
						if(reverseLS[j-1] >= j){
							continue;
						}

						int right = reverseText.length
								- (reverseText.length / 2 - reverseLS[j - 1]) - 1
								- currentText.length() % 2;
						int left = reverseText.length
								- (reverseText.length / 2 + j + reverseLP[j] - 1)
								- 1 - currentText.length() % 2;

						if (currentTreeDepth == 0
								|| !(left == 0 && currentLeft != 0)
								&& !(right == currentText.length() - 1 && (currentLeft + currentText
										.length()) != currentZeroDepth)) {
							left += currentLeft;
							right += currentLeft;

							repetition_list.addRepetition(j, left, right);
						}

					}				
				}
				if (currentTreeDepth == 0) {
					localStack.push(new StackRecord(u, currentLeft, currentTreeDepth+1, currentText.length()));
				} else {
					localStack.push(new StackRecord(u, currentLeft, currentTreeDepth+1, currentZeroDepth));
				}

				if (currentTreeDepth == 0) {
					localStack.push(new StackRecord(v, currentLeft+u.length(), currentTreeDepth+1, currentText.length()));
				} else {
					localStack.push(new StackRecord(v, currentLeft+u.length(), currentTreeDepth+1, currentZeroDepth));
				}
			}
		}
	}

*/
	public void extendMaximalRepetition(String id, String text, RepetitionList repetition_list) {
		int maxRep = 0;
		int maxRepi = -1;


		for(int i = 0; i < repetition_list.get_Repetition_ListSize(); i++){
			int len = repetition_list.get_MaximalRepetition_length(i);
			if(len> maxRep){
				maxRepi = i;
				maxRep = len;
			}
		}

		if(maxRepi == -1) {
			return;
		}

		// repetition_list.print_repetition(text);

		if (maxRep >= 8){
			extendSeed(id, text, repetition_list, maxRepi);
		}
		//	    repetition_list.print_repetition(text);
	}

	public void extendSeed(String id, String text, RepetitionList repetition_list, int index) {
		int left = repetition_list.get_repetition_left_index(index);
		int right = repetition_list.get_repetition_right_index(index);

		// System.out.println(id);

		// System.out.println(text + " " + left+" "+right);
		// leftｓ
		int l = extendLeft(text, repetition_list, index, left);
		int r = extendRight(text, repetition_list, index, right);
		// System.out.println(l+" "+r);

		if(l == left && r == right){
			return;
		}

		int unit = repetition_list.get_repetition_unit_size(index);

		repetition_list.insert(unit, l, r);

	}

	private int extendLeft(String text, RepetitionList repetition_list, int index, int left) {
		int start = left -2;
		int unit = repetition_list.get_repetition_unit_size(index);
		//int startU = (left -start)  % unit;
		if(start <= 0){
			return start + 2;
		}	
		while(true){
			int mis = extendLeft1Mismatch(text, left, start, unit);
			int ins = extendLeft1Ins(text, left, start, unit);
			int del = extendLeft1Del(text, left, start, unit);
			//  System.out.println(mis+"mis "+ins+" "+del);
			if (mis <= ins ) {
				if (mis <= del){ 
					if (start - mis >= 3 && start - mis + 1>= unit ){
						start = mis -2 ;
						continue;
					}else {
						return start + 2; 
					}
				}else {
					if (start-del >= 3 && start -del + 1 >= unit ){
						start = del - 2 ;
						continue;
					}else {
						return start + 2; 
					}
				}
			}else if (ins <= del){
				if (start -ins >= 3 && start -ins + 1 >= unit ){
					start =  ins -2;
					continue;
				}else {
					return start + 2; 
				}
			}else {
				if (start -del>= 3 && start -del + 1 >= unit ){
					start = del -2;
					continue;
				}else {
					return start + 2; 
				}
			}
		}
	}

	private  int extendLeft1Mismatch(String text, int left, int start, int unit) {
		int i;
		for(i= start; i>=0; i--){
			int  Uindex = (i-left)% unit +unit;
			if (text.charAt(i) != text.charAt(left+Uindex)){
				break;
			}
		}
		return ++i;
	}

	private  int extendLeft1Ins(String text, int left, int start, int unit) {
		int i;
		for(i= start; i>=0; i--){
			int  Uindex = (i -left +1)  % unit +unit;
			if (text.charAt(i) != text.charAt(left+Uindex)){
				break;
			}
		}
		return ++i ;
	}

	private  int extendLeft1Del(String text, int left, int start, int unit) {
		int i;
		for(i= start+1; i>=0; i--){
			int  Uindex = (i-left -1)  % unit +unit;
			if (text.charAt(i) != text.charAt(left+Uindex)){
				break;
			}
		}
		return ++i;
	}



	private int  extendRight(String text, RepetitionList repetition_list, int index, int right) {
		int start = right+ 2;
		int unit = repetition_list.get_repetition_unit_size(index);
		//int startU = (right -start)  % unit;
		if(start >= text.length()){
			return start - 2;
		}
		while(true){
			int mis = extendRight1Mismatch(text, right, start, unit);
			int ins = extendRight1Ins(text, right, start, unit);
			int del = extendRight1Del(text, right, start, unit);
			//System.out.println(mis+" "+ins+" "+del);
			if (mis >= ins ) {
				if (mis >= del){

					if (mis -start  >= 3 && mis - start + 1 >= unit ){
						start = mis +2;
						continue;
					}else {
						return start -2 ; 
					}
				}else {
					if ( del-start >= 4 && del -start + 1  >= unit ){
						start = del + 2 ;
						continue;
					}else {
						return start-2; 
					}
				}
			}else if (ins >= del){
				if (ins - start >= 3 && ins - start + 1 >= unit ){
					start = ins + 2;
					continue;
				}else {
					return start -2; 
				}
			}else {
				if ( del -start >= 3 && del -start + 1 >= unit ){
					start = del +2;
					continue;
				}else {
					return start-2; 
				}
			}
		}
	}


	private  int extendRight1Mismatch(String text, int right, int start, int unit) {
		int i;
		for(i= start; i<text.length(); i++){
			int  Uindex = (i-right)% unit -unit;
			if (text.charAt(i) != text.charAt(right+Uindex)){
				break;
			}
		}
		return --i;
	}

	private  int extendRight1Ins(String text, int right, int start, int unit) {
		int i;
		for(i= start; i<text.length(); i++){
			int  Uindex = (i -right -1)  % unit -unit ;
			if (text.charAt(i) != text.charAt(right+Uindex)){
				break;
			}
		}
		return --i ;
	}

	private  int extendRight1Del(String text, int right, int start, int unit) {
		int i;
		for(i= start+1; i<text.length(); i++){
			int  Uindex = (i-right +1)  % unit -unit;
			if (text.charAt(i) != text.charAt(right+Uindex)){
				break;
			}
		}
		return --i;
	}






	private int[] calc_lppattern2(String pattern) {

		int[] lppattern = new int[pattern.length()];

		lppattern[0] = pattern.length();

		for (int i = 0; i < pattern.length(); i++) {
			int j = 0;
			while (pattern.charAt(j) == pattern.charAt(i + j)) {
				j++;
				if (i + j == pattern.length()) {
					break;
				}
			}
			lppattern[i] = j;
		}

		return lppattern;
	}

	private int[] calc_lptext2(String text, String pattern) {
		int[] lptext = new int[text.length()];

		lptext[0] = text.length();

		for (int i = 0; i < text.length(); i++) {
			int j = 0;
			while (text.charAt(i + j) == pattern.charAt(j)) {
				j++;
				if (i + j == text.length() || j == pattern.length()) {
					break;
				}
			}
			lptext[i] = j;
		}

		return lptext;
	}

	public String generateRandomNumbers(int targetLen) {
		Random rnd = new Random();
		char[] target = new char[targetLen];
		// Generate an array in which elements are selected from {1,2,3,4} at
		// random.
		for (int i = 0; i < targetLen; i++) {
			int a = rnd.nextInt(5);
			switch (a) {
			case 0:
				target[i] = 'A';
				break;
			case 1:
				target[i] = 'C';
				break;
			case 2:
				target[i] = 'G';
				break;
			case 3:
				target[i] = 'T';
				break;
			case 4:
				target[i] = 'N';
				break;
			}
		}

		return String.valueOf(target);
	}
}
