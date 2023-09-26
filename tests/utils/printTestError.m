function printTestError(filename, exception)
    fprintf('\n********************************************\n')
    fprintf('****** Dynare test error encountered *******\n')
    fprintf('********************************************\n')
    fprintf('  While running file: %s\n', filename);
    if isoctave
        fprintf('             message: %s\n', exception.message)
        if (isfield(exception, 'stack'))
            fprintf('\n  Stacktrace:\n')
            for i=1:size(exception.stack, 1)
                fprintf('  %s(%s):%d.%d\n', exception.stack(i).file, exception.stack(i).name, ...
                        exception.stack(i).line, exception.stack(i).column);
            end
        end
    else
        fprintf('\n')
        disp(getReport(exception))
    end
    fprintf('********************************************\n\n\n')
end
