package domain.io;

import java.io.File;
import java.io.FileFilter;

public class Filter implements FileFilter {
	private final String extention; 
	public Filter(String ext){
		extention = ext;
	}
	public boolean accept(File pathname) {
		return pathname.getName().endsWith(extention);  
	}

}
