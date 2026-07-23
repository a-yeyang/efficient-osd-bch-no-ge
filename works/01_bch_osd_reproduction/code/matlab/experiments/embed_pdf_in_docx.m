function info = embed_pdf_in_docx(options)
%EMBED_PDF_IN_DOCX Add a PDF OLE object (when Word supports it) or hyperlink.
%
% This optional Windows/Word-COM entry is a safe MATLAB counterpart of the
% Python DOCX ZIP post-processor.  It never overwrites its source document by
% default: it writes assets/generated/*_with_reference_paper.docx.  Word's
% installed PDF OLE handler is environment-specific; failure to create an OLE
% icon transparently falls back to a working hyperlink and reports that fact.

if nargin < 1, options=struct(); end
paths=work01_setup();
defaultInput=fullfile(paths.assets_root,'generated','Work01_MATLAB_Reproduction_Report.docx');
inputPath=local_opt(options,'input_path',defaultInput);
defaultOut=fullfile(paths.assets_root,'generated','Work01_MATLAB_Reproduction_Report_with_reference_paper.docx');
outPath=local_opt(options,'output_path',defaultOut);
visible=local_opt(options,'visible',false); overwrite=local_opt(options,'overwrite',false);
pdfPath=local_opt(options,'pdf_path',fullfile(paths.docs_dir,'reference_paper', ...
    'Efficient_Ordered_Statistics_Decoding_of_BCH_Codes_Without_Gaussian_Elimination.pdf'));
if ~ispc, error('work01:DOCX:WindowsOnly','Word COM automation is available only on Windows.'); end
if ~isfile(inputPath), error('work01:DOCX:MissingInput','Input DOCX not found: %s',inputPath); end
if ~isfile(pdfPath), error('work01:DOCX:MissingPDF','Reference PDF not found: %s',pdfPath); end
if ~isfolder(fileparts(outPath)), mkdir(fileparts(outPath)); end
if isfile(outPath) && ~overwrite, outPath=local_unique_path(outPath); end

word=[]; doc=[]; embedded=false; fallback=false;
try
    word=actxserver('Word.Application'); word.Visible=logical(visible);
    doc=word.Documents.Open(inputPath);
    range=doc.Content; range.Collapse(0); % wdCollapseEnd
    range.InsertAfter(sprintf('\r\nReference paper PDF: '));
    range.Collapse(0);
    try
        % Empty ClassType asks Word to select the locally registered PDF OLE server.
        doc.InlineShapes.AddOLEObject('',pdfPath,false,true,'',0, ...
            'Reference paper PDF',range);
        embedded=true;
    catch
        fallback=true;
        try
            doc.Hyperlinks.Add(range,pdfPath,[],[],'Open reference-paper PDF');
        catch
            range.InsertAfter(pdfPath);
        end
    end
    doc.SaveAs2(outPath,16); doc.Close(false); doc=[]; word.Quit; delete(word); word=[];
    info=struct('supported',true,'output_path',outPath,'pdf_path',pdfPath, ...
        'ole_embedded',embedded,'hyperlink_fallback',fallback,'used_word_com',true);
catch ex
    local_cleanup_word(doc,word); rethrow(ex);
end
end

function local_cleanup_word(doc,word)
try, if ~isempty(doc), doc.Close(false); end, catch, end
try, if ~isempty(word), word.Quit; delete(word); end, catch, end
end
function pathOut=local_unique_path(pathIn)
[folder,name,ext]=fileparts(pathIn); pathOut=pathIn; suffix=1;
while isfile(pathOut), pathOut=fullfile(folder,sprintf('%s_%d%s',name,suffix,ext)); suffix=suffix+1; end
end
function value=local_opt(s,name,default)
if isfield(s,name), value=s.(name); else, value=default; end
end
