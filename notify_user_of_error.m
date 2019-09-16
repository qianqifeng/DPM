function              notify_user_of_error(err)

    clear_waitbars();

    switch(err.identifier)
        case 'MATLAB:unassignedOutputs'
            fprintf('DPM:Model function error \n \t Make sure all the output arguments are set in the \n\t model function (X, C, I, out).\n')
        case 'MATLAB:class:SetProhibited'
            fprintf('DPM:Model function error \n \t Make sure all the output arguments are set in the \n\t model function (X, C, I, out).\n')
            fprintf('  or \n \t Make sure the each of the output states X{1}, X{2},..., X{Nx} have \n\t the same dimensions as input states inp.X{1}, inp.X{2},..., inp.X{Nx}.\n')
            fprintf('  or \n \t Make sure the each of the outputs C{1} and I have \n\t the same dimensions as the input states inp.X{1}, inp.X{2},..., inp.X{Nx}.\n')
        case 'MATLAB:cellRefFromNonCell'
            fprintf('DPM:Model function error \n \t Make sure all the output C is a 1x1 cell array and \n\t that the output X is a cell array with as many elements \n\t as state variables.\n')
        case 'MATLAB:dimagree'
            fprintf('DPM:Model function error \n \t Make sure the each of the output states X{1}, X{2},..., X{Nx} have \n\t the same dimensions as input states inp.X{1}, inp.X{2},..., inp.X{Nx}.\n')
        case 'DPM:Internal'
            %fprintf(['DPM:Error \n \t' err.message '\n'])
            fprintf('DPM:Error \n \t Check the model function.\n')
            rethrow(err)
        otherwise
%             err.identifier
            fprintf('DPM:Error \n \t Check the model function.\n')
            rethrow(err)
    end
    function clear_waitbars()
        
    % clear all waitbars
    set(0,'ShowHiddenHandles','on')
    handles = get(0,'Children');
    for i=1:length(handles)
        if strcmp(get(handles(i),'name'),'DPM:Waitbar')
            delete(handles(i))
        end
    end
    set(0,'ShowHiddenHandles','off')        
