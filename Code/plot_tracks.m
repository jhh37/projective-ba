function h = plot_tracks(M, W, varargin)
% PLOT_TRACKS
% For measurement matrix coming from 2D point tracks, plot the tracks.

opts = au_opts( ...
  'show_timeline=0', ...
  'time_interval=0.3', ...
  'x_min=0', ...
  'x_max=800', ...
  'y_min=0', ...
  'y_max=800', ...
  'negative_measurements=0', ...
  'negative_negative=0', ...
  'save_images=0', ...
  varargin{:});

if opts.negative_negative, M = -M;
end

X = M(1:2:end, :);
Y = M(2:2:end, :);
W2 = 1 - W;
W(W==0) = nan;
W2(W2==0) = nan;

if issparse(Y)
  Y = full(Y);
  Y(Y==0) = nan;
end

if ~opts.show_timeline
  % If not showing a timeline, just plot the measurement matrix.
  h = plot(X, Y);
else
  % Otherwise, show the time evolution of the tracking points.
  T = size(X, 1);
  % h = figure;
  while 1
    for t = 1 : T
      h = plot(W2(2 * t - 1, :) .* X(t, :), W2(2 * t, :) .* Y(t, :), 'LineStyle', 'none', 'Marker', 'o', 'MarkerSize', 3, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');      
      hold on
      h = plot(W(2 * t - 1, :) .* X(t, :), W(2 * t, :) .* Y(t, :), 'LineStyle', 'none', 'Marker', 'o', 'MarkerSize', 3, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');      
      hold off
      if ~opts.negative_measurements
        % axis([opts.x_min opts.x_max opts.y_min opts.y_max]);
        axis([-800 800 -800 800]);
      else
        axis([-opts.x_max opts.x_min -opts.y_max opts.y_min]);
      end
      
      axis off
      
      if opts.save_images
        print('-djpeg', sprintf('%04d.jpg', t));
      end
      
      pause(opts.time_interval);
      clf(h);      
    end
    
    % Turn off the figure-saving option.
    opts.save_images = 0;
  end
  
  % if 0
  %   h(size(M,2))=0;
  %   holdstatus = ishold;
  %   for k=1:size(M,2)
  %     i = find(M(:,k));
  %     col = M(i,k);
  %     x = col(1:2:end);
  %     y = col(2:2:end);
  %
  %     h(k) = plot(x, y, varargin{:});
  %     hold on
  %   end
  %   if ~holdstatus, hold off; end
  
  
end
