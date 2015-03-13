package main;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.StringTokenizer;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.TimeUnit;

import domain.CustomThreadPoolExecutor;
import domain.Filter;
import domain.GenomicLocation;
import domain.KmerRepeatUnitPair;
import domain.kmerInformator;
import domain.multiThreading.processingThread;

public class MainDriver {

	public static void main(String[] args) throws IOException, InterruptedException {
		File directory = new File("/media/gabdank/Disk3/mySTR/MY2");
		File[] files = directory.listFiles(new Filter("fastq"));  
		for (int index = 0; index < files.length; index++) {  
			String filePath = files[index].toString();
			String fileName = "ll";
			StringTokenizer t = new StringTokenizer(filePath,"/");
			while (t.hasMoreTokens()){
				fileName = t.nextToken();
			}
			if (fileName.contains("_R1_")){
				System.out.println(fileName);
				System.out.println(fileName.replace("_R1_", "_R2_"));
				System.out.println();
			}

		} 

		int numberOfThreads = files.length/2;
		System.out.println("Number of threads is: "+numberOfThreads);
		System.out.println("Number of parralel jobs is: "+6);//args[0]
		int parallel = 6;//Integer.parseInt(args[0]);

		long time1 = System.currentTimeMillis();

		kmerInformator testInformator = new kmerInformator(12, "/media/gabdank/Disk3/mySTR/chromosomes.list", "/media/gabdank/Disk3/mySTR/WS243.indexed.bed", 13);

		long time2 = System.currentTimeMillis();
		Integer threadCounter = 0;
		BlockingQueue<Runnable> blockingQueue = new ArrayBlockingQueue<Runnable>(200);
		CustomThreadPoolExecutor executor = new CustomThreadPoolExecutor(parallel,
				parallel, 1, TimeUnit.MILLISECONDS, blockingQueue);

		// Let start all core threads initially
		executor.prestartAllCoreThreads();
		int threadsCounter = 0;
		for (int index = 0; index < files.length; index++)  
		{  
			String filePath = files[index].toString();
			String fileName = "ll";
			StringTokenizer t = new StringTokenizer(filePath,"/");
			while (t.hasMoreTokens()){
				fileName = t.nextToken();
			}
			if (fileName.contains("_R1_")){
				String file1 = fileName;
				String file2 = fileName.replace("_R1_", "_R2_");
				executor.execute(new processingThread(threadsCounter,directory+"/"+file1,directory+"/"+file2, testInformator,directory+"/output_") );
				threadsCounter++;		       

			}
		}  


		executor.shutdown();
		while(! executor.isTerminated()){
			Thread.sleep(100000);
		}

		long time3 = System.currentTimeMillis();

		System.out.println("The genomic locations data structure loading took : "+(time2-time1));
		System.out.println("The multithreaded part took : " +(time3-time2));
		System.out.println("The run in total took : "+(time3-time1));




	}

}
