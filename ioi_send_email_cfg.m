function sendmail1 = ioi_send_email_cfg
% Graphical interface configuration that sends an automatic e-mail notifier.
%_______________________________________________________________________________
% Copyright (C) 2012 LIOM Laboratoire d'Imagerie Optique et Mol�culaire
%                    �cole Polytechnique de Montr�al
%_______________________________________________________________________________

% ------------------------------------------------------------------------------
% Recipient
% ------------------------------------------------------------------------------
recipient           = cfg_entry;
recipient.tag       = 'recipient';
recipient.name      = 'Recipient';
recipient.val       = {''};
recipient.strtype   = 's';
recipient.num       = [1 Inf];
recipient.help      = {'User to receive mail.'};
% ------------------------------------------------------------------------------
% Subject
% ------------------------------------------------------------------------------
subject             = cfg_entry;
subject.tag         = 'subject';
subject.name        = 'Subject';
subject.val         = {['[ioi11] [%DATE%] On behalf of ' spm('Ver')]};
subject.strtype     = 's';
subject.num         = [1 Inf];
subject.help        = {'The subject line of the message. %DATE% will be replaced by a string containing the time and date when the email is sent.'};
% ------------------------------------------------------------------------------
% Message
% ------------------------------------------------------------------------------
message             = cfg_entry;
message.tag         = 'message';
message.name        = 'Message';
message.val         = {['Hello from ioi11!' sprintf('\n') 'Job finished on %DATE%' ]};
message.strtype     = 'e';
message.num         = [1 Inf];
message.help        = {'A string containing the message to send. %DATE% will be replaced by a string containing the time and date when the email is sent.'};
% ------------------------------------------------------------------------------
% Attachments
% ------------------------------------------------------------------------------
attachments         = cfg_files;
attachments.tag     = 'attachments';
attachments.name    = 'Attachments';
attachments.val{1}  = {};
attachments.filter  = '.*';
attachments.ufilter = '.*';
attachments.num     = [0 Inf];
attachments.help    = {'List of files to attach and send along with the message.'};
% ------------------------------------------------------------------------------
% SMTP Server
% ------------------------------------------------------------------------------
smtp                = cfg_entry;
smtp.tag            = 'smtp';
smtp.name           = 'SMTP Server';
smtp.strtype        = 's';
smtp.val            = {'smtp.gmail.com'};
smtp.num            = [1 Inf];
smtp.help           = {'Your SMTP server. If not specified, look for sendmail help.'};
% ------------------------------------------------------------------------------
% E-mail
% ------------------------------------------------------------------------------
email               = cfg_entry;
email.tag           = 'email';
email.name          = 'E-mail';
email.strtype       = 's';
email.val           = {'ioi11notifier@gmail.com'};
email.num           = [1 Inf];
email.help          = {'The e-mail address of the notifier. Look in sendmail help how to store it. Default: ioi11notifier@gmail.com'};
% ------------------------------------------------------------------------------
% Password
% ------------------------------------------------------------------------------
password            = cfg_entry;
password.tag        = 'password';
password.name       = 'password';
password.strtype    = 's';
password.val        = {'epoxy111'};
password.num        = [1 Inf];
password.help       = {'Your e-mail password. Look in sendmail help how to store it.'};
% ------------------------------------------------------------------------------
% Zip attachments
% ------------------------------------------------------------------------------
zip                 = cfg_menu;
zip.tag             = 'zip';
zip.name            = 'Zip attachments';
zip.val             = {'No'};
zip.labels          = {'Yes' 'No'}';
zip.values          = {'Yes' 'No'}';
zip.help            = {'Zip attachments before being sent along with the message.'};
% ------------------------------------------------------------------------------
% Parameters
% ------------------------------------------------------------------------------
params              = cfg_branch;
params.tag          = 'params';
params.name         = 'Parameters';
params.val          = { smtp email password zip };
params.help         = {'Preferences for your e-mail server (Internet SMTP server) and your e-mail address. If you encounter any error, identify the outgoing mail server for your electronic mail application, which is usually listed in the application''s preferences, or, consult your e-mail system administrator, and update the parameters. Note that this function now supports e-mail servers that require authentication, like gmail.'};
% ------------------------------------------------------------------------------
% Sendmail
% ------------------------------------------------------------------------------
sendmail1           = cfg_exbranch;
sendmail1.tag       = 'sendmail1';
sendmail1.name      = 'Send e-mail notifier';
sendmail1.val       = { recipient subject message attachments params };
sendmail1.prog      = @ioi_sendmail;
sendmail1.help      = {'Send a mail message (attachments optionals) to an address.'};
end % ioi_send_email_cfg
%_______________________________________________________________________________

%_______________________________________________________________________________
function ioi_sendmail(job)
% Sets up SMTP server and sends e-mail.
try
    % Set up the preferences properly:
    setpref('Internet','E_mail',job.params.email);
    setpref('Internet','SMTP_Server',job.params.smtp);
    setpref('Internet','SMTP_Username',job.params.email);
    setpref('Internet','SMTP_Password',job.params.password);
    props = java.lang.System.getProperties;
    props.setProperty('mail.smtp.auth','true');
    props.setProperty('mail.smtp.socketFactory.class', 'javax.net.ssl.SSLSocketFactory');
    props.setProperty('mail.smtp.socketFactory.port','465');
    
    % Format the subject line & message
    subj = strrep(job.subject,'%DATE%',spm('time'));
    mesg = strrep(job.message,'%DATE%',spm('time'));
    mesg = [mesg sprintf('\n') ...
    '%_______________________________________________________________________________' sprintf('\n') ...
    '% Copyright (C) 2012 LIOM Laboratoire d''Imagerie Optique et Mol�culaire' sprintf('\n') ...
    '%                    �cole Polytechnique de Montr�al' sprintf('\n') ...
    '%_______________________________________________________________________________'...
    sprintf('\n --  \n Statistical Parametric Mapping')];
    
    % Prepare attachments
    if ~isempty(job.attachments)
        if strcmpi(job.params.zip,'Yes')
            zipfile = fullfile(tempdir,'ioi_sendmail.zip');
            zip(zipfile,job.attachments);
            job.attachments = {zipfile};
        end
        sendmail(job.recipient,subj,mesg,job.attachments);
    else
        sendmail(job.recipient,subj,mesg);
    end
    fprintf('Sending mail %60s...\n', spm('time'));
catch
    %- not an error to prevent an analysis to crash because of just that...
    fprintf('Sendmail failed...\n');
end
end % ioi_sendmail

% EOF
