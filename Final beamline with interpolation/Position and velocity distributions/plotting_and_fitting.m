for i = 1:n_out_of_beam_source
    v_final(i,:,:) = v_out_of_beam_source{i};
    x_final(i,:,:) = x_out_of_beam_source{i};
end


%Make histograms of velocities and positions:
figure
hist(v_final(:,1,n_end))
title('v_x after exiting beambox')

figure
hist(v_final(:,2,n_end))
title('v_y after exiting beambox')

figure
hist(v_final(:,3,n_end))
title('v_z after exiting beambox')

%Plot positions after beam box
figure
n = hist3(x_final(:,[1 2],n_end),'Nbins',[200 200]);
pcolor(n)

%% save the final positions and velocities in a file:
save('velocities_after_beamsource_1e6.mat','v_final')
save('positions_after_beam_source_1e6.mat','x_final')

