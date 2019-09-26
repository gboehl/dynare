function printMakeCheckOctaveErrMsg(modfilename, err)
    printf("\n");
    printf("********************************************\n");
    printf("*** DYNARE-TEST-OCTAVE ERROR ENCOUNTERED ***\n");
    printf("********************************************\n");
    printf("  WHILE RUNNING MODFILE: %s\n", modfilename);
    printf("                    MSG: %s\n", err.message);
    if (isfield(err, 'stack'))
        printf("\n  STACKTRACE:\n");
        for i=1:size(err.stack, 1)
            printf("  %s(%s):%d.%d\n", err.stack(i).file, err.stack(i).name,
                   err.stack(i).line, err.stack(i).column);
        endfor
    end
    printf("********************************************\n\n\n");
end
