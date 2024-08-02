laraImport("./lara/StressTest");

function Tester(repetitions = 3) {
    console.log("Calling StressTest " + repetitions + " times");
    for (let i = 0; i < repetitions; i++) {
        const index = i + 1;
        console.log("StressTest " + index + "/" + repetitions);
        const statsFilename = "stats_" + index + ".json";

        StressTest("src", statsFilename);
        console.log("Results saved to " + statsFilename);
    }
}
