package com.manutils.vcfBedIntersection;

import IntervalSearch.*;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class FindIntersection {

//	static String vcfPath = "C:\\Users\\pc\\Desktop\\VIkasCOde\\filesandscript\\SNVs_Test.vcf";
//	static String bedPath = "C:\\Users\\pc\\Desktop\\VIkasCOde\\filesandscript\\TargetedRegion.bed";
//	static String outFilepath = "C:\\Users\\pc\\Desktop\\VIkasCOde\\filesandscript\\testOtDebug.txt";
	static Map<String, IntervalTree<Long>> hmIt = new HashMap<String, IntervalTree<Long>>();

	public static void main(String args[])throws IOException{

		String vcfPath = args[0];
		String bedPath = args[1];
		String outFilepath = args[2];
		
		BufferedReader br = new BufferedReader(new FileReader(bedPath));
		PrintWriter pw = new PrintWriter(new File(outFilepath));
		String record = ""; 
		String[] splitedRecordArr;
		IntervalTree<Long> it = new IntervalTree<Long>();
		Long chrCounter = 0L;
		String previousChromosome="";
		int recordsCounter = 0;
		while((record = br.readLine())!=null){
			if(record.startsWith("chr")){

				splitedRecordArr = record.split("\\t",4);
				Long start = Long.parseLong(splitedRecordArr[1]);
				Long stop = Long.parseLong(splitedRecordArr[2]);

				if(recordsCounter == 0){
					recordsCounter++;
					it.addInterval(start, stop, chrCounter);
					previousChromosome = splitedRecordArr[0];
					continue;
				}

				if((splitedRecordArr[0].equals(previousChromosome))){
					chrCounter++;
					it.addInterval(start, stop, chrCounter);

					hmIt.put(previousChromosome, it);
				}else if(!(splitedRecordArr[0].equals(previousChromosome))){
					hmIt.put(previousChromosome, it);
					chrCounter=0L;
					it = new IntervalTree<Long>();
					it.addInterval(start, stop, chrCounter);
					hmIt.put(splitedRecordArr[0], it);
				}

				previousChromosome = splitedRecordArr[0];
			}
		}
		br.close();

		BufferedReader brVcf = new BufferedReader(new FileReader(vcfPath));

		while((record = brVcf.readLine()) != null){
			if(record.startsWith("chr")){
				splitedRecordArr = record.split("\\t",3);
				String chrm = splitedRecordArr[0];
				Long pos = Long.parseLong(splitedRecordArr[1]);
				if(hmIt.containsKey(chrm)){
					List<Long> l = hmIt.get(chrm).get(pos);

					if(!l.isEmpty()){
						pw.append(record+"\n");
					}
				}	
			}
		}
		brVcf.close();
		pw.close();
	}

}
