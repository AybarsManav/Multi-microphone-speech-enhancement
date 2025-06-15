function player = play_scaled_sound(sound, fs)
    sound = normalize(sound, "range", [-1, 1]);
    player = audioplayer(sound, fs);
    play(player);
end

