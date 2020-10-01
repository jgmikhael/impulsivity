% This code simulates the effect of dopamine on behavior in intertemporal
% choice tasks both with (van Gaalen et al., 2006) and without (Cardinal et
% al., 2000; Tanno et al., 2014) option-specific post-reward delays.
% Finally, it reproduces the saline conditions in Fig. 3B,C of 
% Cardinal et al. (2000) within the same figure.
% Written 14Aug20 by JGM.

%-------------------------------------------------------------------------%

% parameters    
alphaT = .15;             	% Weber fraction for time
dag = [0 .5 1.5 5];      	% DA agonist in arbitrary units
d = 1+dag;              	% total DA level
beta = 25;                 	% inverse temperature for choice rule
alphaR = .4;              	% Weber fraction for rewards
epsR = .5;                  % signal-independent noise for rewards
epsT = 1;                   % signal-independent noise for time
POSTr = 50;               	% arbitrary POST if no contingency
g = [1; 8];              	% gain on temporal SD; [high; low] precision
g = [g; g(1)];            	% [high low high] for figures 1, 2, 3
eta = 8;                   	% POST relative pacemaker period

%-------------------------------------------------------------------------%

% pre-reward delay, for [small; large]
PREx = 0:1:60;                          % PRE for larger/later
PRE = [0*PREx; PREx];                   % likelihood means
PRE0 = mean(PRE);                       % prior mean
PREl0 = 1./(1+std(PRE).^2);             % prior precision
tt = .75*(-5:PREx(end));                % time domain for illustration

% rewards, for [small large]
r = [1 4];                              % likelihood means
r0 = mean(r);                           % prior mean
rl0 = 1./(1+std(r)^2);                  % prior precision
rr = 0:.01:1.5*max(r);                  % reward domain for illustration

figName = {'main1','main2','VV06b','baseline','saline'};
figInd = 0;                             % initialize figure numbers
C = [0 .3 .6 .8]'*[1 1 1];              % color scheme
sampl = 1:10:(max(PRE(:))+1);           % marker position
msize = 120;                            % marker size
mrkr = {'o','s','s','^'};               % marker type

figure(101)
for expt = 1:3                          % representing figures 1,2,3

    % POST
    if expt < 3                         % no option-POST contingency
        POST = 0*PRE + POSTr/eta;       % likelihood mean
    elseif expt == 3                    % with option-POST contingencies
        POST = [100+0*PREx; 100-PREx]/eta;
    end
    POST0 = mean(POST);                 % prior mean
    POSTl0 = 1./(1+std(POST).^2);       % prior precision
    
    p = zeros(length(d),length(PRE));	% initialize p(LL)
    for q = 1:length(d)   

        % effect of DA on likelihood standard deviations
        rs = (epsR+alphaR*r)./d(q);                         % reward
        PREs = (epsT+alphaT*PRE).*g(expt)./d(q);           	% PRE
        POSTs = eta*(epsT+alphaT*POST).*g(expt)./d(q);     	% POST
        
        % reward
        rl = 1./rs.^2;                              % likelihood precisions
        rlh = rl+rl0;                            	% posterior precisions
        rh = (rl.*r+rl0.*r0)./rlh;                  % posterior means
        
        % PRE
        PREl = 1./PREs.^2;                          % likelihood precisions
        PRElh = PREl+PREl0;                      	% posterior precisions
        PREh = (PREl.*PRE+PREl0.*PRE0)./PRElh;    	% posterior means
        
        % POST
        POSTl = 1./POSTs.^2;                        % likelihood precisions
        POSTlh = POSTl+POSTl0;                     	% posterior precisions
        POSTh = (POSTl.*POST+POSTl0.*POST0)./POSTlh;% posterior means
        
        % reward rates
        RS = rh(1)./(PREh(1,:)+POSTh(1,:));       	% small
        RL = rh(2)./(PREh(2,:)+POSTh(2,:));        	% large
        
        % p(selecting large reward)
        p(q,:) = 1./(1+exp(-beta*(RL-RS)));        	% softmax
 
        if q == 1 || q == length(d)                 % lowest, highest DA
            figure(101)
            subplot(6,2,1+4*(expt-1))
            x = normpdf(rr,rh',1./sqrt(rlh'));
            plot(rr,x,'Color',C(q,:))
            hold on
            plot(rh(1)*[1 1],[0 10],'Color',C(q,:),'LineWidth',2)
            hold on
            plot(rh(2)*[1 1],[0 10],'Color',C(q,:),'LineWidth',2)
            hold on
            xlabel('Reward (Numerator)')
            ylabel('p(reward)')
            box off
            % ylim([0 max(x(:))])
            ylim([0 1])
            
            subplot(6,2,3+4*(expt-1))          % illustrate for middle time
            tmid = PREh(:,(end-1)/2);
            SDmid = 1./sqrt(PRElh(:,(end-1)/2));
            x = normpdf(tt,tmid,SDmid);
            plot(tt,x,'Color',C(q,:))
            hold on
            plot(tmid(1)*[1 1],[0 10],'Color',C(q,:),'LineWidth',2)
            hold on
            plot(tmid(2)*[1 1],[0 10],'Color',C(q,:),'LineWidth',2)
            hold on
            xlabel('Time (Denominator)')
            ylabel('p(time)')
            box off
            % ylim([0 max(x(:))])
            ylim([0 .25])
        end
    end
    baseline(expt,:) = p(1,:);                  % baseline (saline) levels
    
    for k = 1:2
    if k == 1; subplot(6,2,4*expt+[-2 0])
        else figInd=figInd+1; figure(figInd)
    end

    for e = 1:length(d)                         % for legend
        DAlevels{e} = [num2str(dag(e)) ' a.u.'];
    end
    
    % plot([0 60],[50 50],'k--'); hold on
    
    for e = 1:size(p,1)
        h(e) = plot(PRE(2,:), 100*p(e,:)','Color',C(e,:));
        hold on
        if expt == 3
            scatter(PRE(2,sampl),100*p(e,sampl),...
                msize,C(e,:),mrkr{e},'filled');
            hold on
        end
    end
    aa1 = legend(h,DAlevels,'Box','Off','FontSize',25);
    aa2 = xlabel('delay (sec)');
    aa3 = ylabel('% choice of the large reinforcer');
    ylim([0 100])
    xticks(0:10:60)
    yticks(0:20:100)
    box off
    
    if expt == 3 && k == 2
        aa1.FontSize = 33;
        aa2.FontSize = 30;
        aa3.FontSize = 30;
        set(gca,'FontSize',28);
    end
    end

    if expt == 2
        figure(4)
        plot(PRE(2,:),100*baseline(1,:),'k')
        hold on
        plot(PRE(2,:),100*baseline(2,:),'Color',.5*[1 1 1])
        ylim([0 100])
        xlabel('delay (sec)')
        ylabel('% choice of the large reinforcer')
        xticks(0:10:60)
        yticks(0:20:100)
        legend('High Temporal Precision','Low Temporal Precision',...
            'Location','Northeast','Box','Off','FontSize',25)
        box off
    end
end

%-------------------------------------------------------------------------%

% reproduce the saline conditions in Cardinal et al. (2000), Fig. 3

house = [81.6; 65.2; 48; 44.5; 38.4];           % houselight
NoCue = [76.2; 54.8; 44.8; 39.5; 28.4];         % no cue
Cue = [82; 56.8; 40; 29; 22];                   % cue

x = [0 10 20 40 60];                            % delay
y = 0:10:100;                                   % percent choice of large

figure(5)
plot(x,Cue,'-ko')
hold on
plot(x,NoCue,'-o','Color',.5+[0 0 0])
legend('Cue','No Cue','Box','Off','FontSize',22)
xticks(x)
yticks(y)
xlabel('Delay to large reinforcer (s)')
ylabel('Percent choice of large reinforcer')
xlim([0 60])
ylim([0 100])
box off

figure(101)