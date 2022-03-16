function [vp, response_mapping, instruct] = get_experimentInfo(run_mode, rootDir)
if run_mode == 1
  fprintf('\nParticipant (three characters, e.g. S01)?\n');
  vp = input('EINGABE: ', 's');
    if length(vp)~=3 
       error ('Wrong input!'); end
    fprintf('\nResponse mapping?\n1: y = rund,  m = spitz\n2: y = spitz, m = rund\n');   
    response_mapping = str2num(input('EINGABE (1 oder 2): ', 's'));    
      if ~ismember(response_mapping, [1, 2])
        error('\nUnknown mapping!'); end
else
      vp = 'tmp';
      response_mapping = 1;
end

switch response_mapping
    case  1
        instruct = imread([rootDir, '/instruction_mapping1.png']);
    case  2
        instruct = imread([rootDir, '/instruction_mapping2.png']);
end
instruct = 255-instruct; instruct(instruct == 255) = 127;
end
