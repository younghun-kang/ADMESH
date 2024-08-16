% Button pushed function: ProgressBarButton
function UpdateProgressBarButton(ProgressBarButton, i,n)
% Change button name to "Processing"
% ProgressBarButton.Text = 'Processing...';
% Put text on top of icon
% ProgressBarButton.IconAlignment = 'Right';
% Create waitbar with same color as button
wbar = permute(repmat(ProgressBarButton.BackgroundColor,15,1,200),[1,3,2]);
% Black frame around waitbar
wbar([1,end],:,:) = 0;
wbar(:,[1,end],:) = 0;
% Load the empty waitbar to the button
% ProgressBarButton.Icon = wbar;
% Loop through something and update waitbar
% n = 10;
% for i = 1:n
    % Update image data (royalblue)
    % if mod(i,10)==0 % update every 10 cycles; improves efficiency
    currentProg = min(round((size(wbar,2)-2)*(i/n)),size(wbar,2)-2);
    RGB = wbar;
    RGB(2:end-1, 2:currentProg+1, 1) = 0.0; % (royalblue)
    RGB(2:end-1, 2:currentProg+1, 2) = 0.4470;
    RGB(2:end-1, 2:currentProg+1, 3) = 0.7410;
    
    ProgressBarButton.Icon = RGB;
    drawnow;
    % Pause to slow down animation
    % pause(.1)
    % end
% end
% remove waitbar
% ProgressBarButton.Icon = '';
% Change button name
% ProgressBarButton.Text = 'Process Data';
end