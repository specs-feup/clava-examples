import StressTest;

aspectdef Tester
	input repetitions = 3 end

	println("Calling StressTest " + repetitions + " times");
	for(var i=0; i<repetitions; i++) {
		var index = i + 1;
		println("StressTest " + index + "/" + repetitions);
		var statsFilename = "stats_" + index + ".json";
		
		call StressTest("src", statsFilename);
		println("Results saved to " + statsFilename);
	}

end
