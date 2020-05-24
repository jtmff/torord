function rewrite_state_vars(InputFile, OutputFile, y)
% By Alfonso Bueno-Orovio
%
% A function that takes a filename for a given CellML InputFile and
% replaces all the initial state variables with the contents of vector y,
% storing the result in OutputFile.
%
% InputFile:  string
% OutputFile: string
% y:          real vector (state variables)

% 1) Define variable names
VarName{1}  = '"v"';
VarName{2}  = '"nai"';
VarName{3}  = '"nass"';
VarName{4}  = '"ki"';
VarName{5}  = '"kss"';
VarName{6}  = '"cai"';
VarName{7}  = '"cass"';
VarName{8}  = '"cansr"';
VarName{9}  = '"cajsr"';
VarName{10} = '"m"';
VarName{11} = '"hp"';
VarName{12} = '"h"';
VarName{13} = '"j"';
VarName{14} = '"jp"';
VarName{15} = '"mL"';
VarName{16} = '"hL"';
VarName{17} = '"hLp"';
VarName{18} = '"a"';
VarName{19} = '"iF"';
VarName{20} = '"iS"';
VarName{21} = '"ap"';
VarName{22} = '"iFp"';
VarName{23} = '"iSp"';
VarName{24} = '"d"';
VarName{25} = '"ff"';
VarName{26} = '"fs"';
VarName{27} = '"fcaf"';
VarName{28} = '"fcas"';
VarName{29} = '"jca"';
VarName{30} = '"nca_ss"';
VarName{31} = '"nca_i"';
VarName{32} = '"ffp"';
VarName{33} = '"fcafp"';
VarName{34} = '"xs1"';
VarName{35} = '"xs2"';
VarName{36} = '"Jrel_np"';
VarName{37} = '"CaMKt"';
VarName{38} = '"C1"';
VarName{39} = '"C2"';
VarName{40} = '"C3"';
VarName{41} = '"O"';
VarName{42} = '"I"';
VarName{43} = '"Jrel_p"';
VarName{44} = '"cli"';
VarName{45} = '"clss"';


% 2) Read whole InputFile into cell array
fid = fopen(InputFile,'r');
file = textscan(fid, '%s', 'Delimiter', '\n', 'CollectOutput', true);
fclose(fid);

% 3) Find positions where changes need to be applied
nVars = length(VarName);
ReplacedStrings = false(nVars,1);
for iline = 1:length(file{1})
    for iVar = 1:nVars
        if (contains(file{1}{iline},VarName{iVar}) && contains(file{1}{iline},'initial_value'))
            disp(iVar)
            oldline = file{1}{iline};
            I = strfind(oldline,'value');
            Idel = I+strfind(oldline(I:end),'"')-1;
            newline = [oldline(1:Idel(1)) num2str(y(iVar),'%.6e') oldline(Idel(2):end)];
            file{1}{iline} = newline; % replace string
            ReplacedStrings(iVar) = true;
        end
    end
end

% 4) Write the modified cell array into the text file
I = find(~ReplacedStrings);
if isempty(I)
    fid = fopen(OutputFile,'w');
    for iline = 1:length(file{1})
        fprintf(fid, '%s\n', char(file{1}{iline}));
    end
    fclose(fid);
else
    fprintf('Output file not written >> unreplaced strings:');
    fprintf(' %d',I);
    fprintf('\n');
end
