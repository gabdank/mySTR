package test;

import static org.junit.Assert.*;

import java.io.IOException;

import org.junit.Test;

import domain.kmerInformator;
import domain.multiThreading.processingThread;

public class ThreadTest {

	@Test
	public void test() throws IOException {
		kmerInformator testInformator = new kmerInformator(12, "/home/gabdank/Documents/STR_Attempt/chromosomes.list", "/home/gabdank/Documents/STR_Attempt/temp.indexed.bed", 5);
		processingThread t = new processingThread(5,"/home/gabdank/workspace/mySTR/fastq/R1.fastq","/home/gabdank/workspace/mySTR/fastq/R2.fastq", testInformator) ;
		t.run();
	}
}
