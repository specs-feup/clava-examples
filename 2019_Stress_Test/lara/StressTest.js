laraImport("lara.Io");
laraImport("lara.Strings");
laraImport("lara.System");
laraImport("lara.Platforms");
laraImport("lara.util.StringSet");
laraImport("clava.Clava");
laraImport("clava.cmake.CMaker");

laraImport("./lara/DynamicCallGraph");

function StressTest(srcFoldername, statsFilename = "stats.json") {
    // To replicate published line diff results, uncomment this line
    //_DCG_USE_COMMA_OPERATOR_ = false;

    setDebug();
    const srcFolder = Io.getPath(
        Clava.getData().getContextFolder(),
        srcFoldername
    );
    const inputFolder = Io.mkdir(Clava.getData().getContextFolder(), "input");
    const outputFolder = Io.mkdir(Clava.getData().getContextFolder(), "output");
    const cmakeFolder = Io.mkdir(
        Clava.getData().getContextFolder(),
        "cmake-build"
    );

    // Arguments for some of the programs
    let args = {};
    const gccCode = Io.getPath(srcFolder, "gcc.c");
    const bzip2Input = Io.getPath(outputFolder, "gcc_for_bzip2.c");
    Io.copyFile(gccCode, bzip2Input);
    const gzipInput = Io.getPath(outputFolder, "gcc_for_gzip.c");
    Io.copyFile(gccCode, gzipInput);

    const oggencInput = Io.getPath(
        inputFolder,
        "pservice-bell_daniel_simion.wav"
    );
    const oggencOutput = Io.getPath(
        outputFolder,
        "pservice-bell_daniel_simion.ogg"
    );
    Io.copyFile(gccCode, gzipInput);

    const gccInput = Io.getPath(inputFolder, "gcc_test.c");
    const gccOutput = Io.getPath(outputFolder, "gcc_test.s");

    args["bzip2"] = [bzip2Input.getAbsolutePath(), "-f"];
    args["gzip"] = [gzipInput.getAbsolutePath(), "-f"];
    args["oggenc"] = [
        oggencInput.getAbsolutePath(),
        "-o=" + oggencOutput.getAbsolutePath(),
        "-Q",
    ];
    args["gcc"] = [
        gccInput.getAbsolutePath(),
        "-o",
        gccOutput.getAbsolutePath(),
    ];
    // Object to collect data
    const stats = {};

    for (const srcFile of Io.getFiles(srcFolder, "*.c")) {
        try {
            console.log("Processing file '" + srcFile + "'");

            const benchName = Io.removeExtension(srcFile.getName());
            //if(srcFile.getName() !== "gcc.c") {
            //	continue;
            //}

            // Add program
            Clava.getProgram().addFileFromPath(srcFile);

            // Rebuild tree
            const parsingStart = System.nanos();
            Clava.rebuild();
            const parsingTime = System.toc(parsingStart);

            // Count Clava nodes before instrumentation
            debug("Counting nodes before instrumentation...");
            const clavaNodes =
                Clava.getProgram().children[0].descendants.length;
            debug("Done, found " + clavaNodes + " nodes");

            // Count lines after parsing but before instrumentation
            const linesBeforeDcg = countNonBlankLines(
                Clava.getProgram().children[0].code
            );

            // Save non instrumented file
            const originalFolder = Io.mkdir(outputFolder, "original");
            Clava.getProgram().children[0].write(
                originalFolder.getAbsolutePath()
            );

            // Graph file
            const graphFile = Io.getPath(outputFolder, benchName + ".dot");

            // Apply dynamic call graph
            const dcgStart = System.nanos();
            DynamicCallGraph(Io.getAbsolutePath(graphFile));
            const dcgTime = System.toc(dcgStart);

            // Count lines after parsing, after instrumentation
            const linesAfterDcg = countNonBlankLines(
                Clava.getProgram().children[0].code
            );
            const linesInserted = linesAfterDcg - linesBeforeDcg;

            // Compile program
            const cmaker = new CMaker("callgraph", false);
            cmaker.addCurrentAst().addLibs("m");

            if (Platforms.isWindows()) {
                cmaker
                    .setGenerator("MinGW Makefiles")
                    .setMakeCommand("mingw32-make");
            }

            const cmakelistsFolder = Io.mkdir(cmakeFolder, benchName);
            const buildFolder = Io.mkdir(cmakelistsFolder, "build");
            const exe = cmaker.build(cmakelistsFolder, buildFolder);

            // Prepare command
            let command = [exe.getAbsolutePath()];
            const extraArgs = args[benchName];
            if (extraArgs !== undefined) {
                command = command.concat(extraArgs);
            }

            // Run program
            debug("Running compiled program: " + command.join(" "));
            const executor = new ProcessExecutor();
            executor.setPrintToConsole(false).execute(command);
            debug("Done running");

            // Save instrumented file
            Clava.getProgram().children[0].write(
                outputFolder.getAbsolutePath()
            );

            // Collect stats
            stats[srcFile.getName()] = collectStats(
                srcFile,
                clavaNodes,
                linesInserted,
                parsingTime,
                dcgTime,
                graphFile
            );

            // Clean program
            Clava.getProgram().removeChildren();

            // Write stats
            Io.writeJson(
                Io.getPath(Clava.getData().getContextFolder(), statsFilename),
                stats
            );
        } catch (e) {
            console.log("Skipping " + srcFile + ", there are problems: " + e);

            // Clean program
            Clava.getProgram().removeChildren();
        }
    }

    console.log("Stats:");
    printlnObject(stats);
}

function collectStats(
    srcFile,
    clavaNodes,
    linesInserted,
    parsingTime,
    dcgTime,
    graphFile
) {
    const stats = {};

    stats["Code Lines"] = Strings.asLines(Io.readFile(srcFile)).length;
    stats["Clava Nodes"] = clavaNodes;
    stats["Inserted Lines"] = linesInserted;
    stats["Parsing"] = parsingTime;
    stats["Dynamic Call Graph"] = dcgTime;

    // Call graph stats
    const callGraph = Io.readFile(graphFile);

    let edges = 0;
    const nodes = new StringSet();

    for (const line of Strings.asLines(callGraph)) {
        if (
            Strings.isEmpty(line) ||
            line.startsWith("}") ||
            line.startsWith("digraph dynamic_call_graph")
        ) {
            continue;
        }

        // Count edge
        edges++;

        // Dot follows this format:
        // NODE1 -> NODE2 [label="<a number>"];
        const parts = line.split("->");
        nodes.add(parts[0].trim());

        const parts2 = parts[1].trim().split(" ");
        nodes.add(parts2[0].trim());
    }

    stats["Graph Nodes"] = nodes.values().length;
    stats["Graph Edges"] = edges;

    return stats;
}

function countNonBlankLines(code) {
    let counter = 0;
    for (const line of Strings.asLines(code)) {
        if (Strings.isEmpty(line)) {
            continue;
        }

        // Do not count cases where {} where inserted
        const strippedLine = line.trim();
        if (strippedLine === "{") {
            continue;
        }

        if (strippedLine === "}") {
            continue;
        }

        counter++;
    }
    return counter;
}
