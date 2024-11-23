// Create time
argList args(argc, argv);

Time runTime
(
    "controlDict",
	args.rootPath(),
	args.caseName(),
	"system",
	"constant",
	!args.optionFound("noFunctionObjects")
);
